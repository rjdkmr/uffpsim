/**
 * @brief This class represents a search engine for molecular fingerprints.
 * @author Rajendra Kumar
 */

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

#include "searchEngine.h"

namespace py = pybind11;

FPSearchEngine::FPSearchEngine(const std::string& filename, std::string mode) {
    _fpStore = new FingerprintStore(filename);
    _molIdMaxLength = _fpStore->_molIdMaxLength;
    _fpSize = _fpStore->_fpSize;
    _CFPSize = _fpStore->_CFPSize;
    _molIdOffset =  _fpStore->_molIdOffset;
    _CFPPopCountIndex = _fpStore->_CFPPopCountIndex;
    _fpEndIndex = _molIdOffset + _fpSize;
    _mode = mode;

    if (mode == "disk") {
        _normal_search = &FPSearchEngine::_normal_search_disk;
        _batch_search = &FPSearchEngine::_batch_search_disk;
        _fpStore->loadDataInMemory(true); // load only cluster fps in memory for disk-based search
        _fpStore->initH5GroupsMappingForPopCountBins(); // initialize the mapping of popcount to corresponding h5 groups for popcount bins and clusters for disk-based search
    } else {
        _normal_search = &FPSearchEngine::_normal_search_memory;
        _batch_search = &FPSearchEngine::_batch_search_memory;
        _fpStore->loadDataInMemory(); // load all fps in memory for memory-based search
    }
    

    _div_lookup_table.assign(static_cast<size_t>(_fpSize * 128 + 2), 0.0f);
    for (int i = 0; i <= _fpSize * 128; i++) {
        _div_lookup_table[i + 1] = 1.0f / static_cast<float>(i + 1);
    }
    py::gil_scoped_acquire acquire;
       
}

FPSearchEngine::~FPSearchEngine() {
    if (_fpStore!= nullptr)
        delete _fpStore;
}

void FPSearchEngine::close() {
    _fpStore->freeMemory();
    _fpStore->close();
    delete _fpStore;
    _fpStore = nullptr;
}

uint64_t* FPSearchEngine::prepareQuery(const std::string& query) {
    uint64_t *queryCFp = new uint64_t[_fpStore->_CFPSize];
    utils::BitStrToCompactFPArray(query, "_", queryCFp, _fpStore->_molIdOffset, _fpStore->_CFPSize);
    return queryCFp;
}

std::vector<utils::dt_inner_clusters_fingerprints_maxscore> FPSearchEngine::filterPopcountBins(uint64_t queryPopcount, float threshold) {
    std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore;
    for (size_t i=0; i < _fpStore->_popCountBins.size(); i++) {
        int popCount = _fpStore->_popCountBins[i];
        int maxPopCountWithQuery = std::max(popCount, (int) queryPopcount);
        float maxScore = (float) std::min(popCount, (int) queryPopcount) / maxPopCountWithQuery;
        if (maxScore >= threshold ) {
            utils::dt_inner_clusters_fingerprints *inner_clusters_fingerprints = _fpStore->_fp_inner_clusters_by_popcount + i;
            utils::dt_inner_clusters_fingerprints_maxscore inner_clusters_fingerprints_maxscore = {*inner_clusters_fingerprints, maxScore, maxPopCountWithQuery};
            filteredPopCountBinsWithMaxScore.push_back(inner_clusters_fingerprints_maxscore);
        }
    }

    // Sort by score
    std::sort(filteredPopCountBinsWithMaxScore.begin(), filteredPopCountBinsWithMaxScore.end(), [](const utils::dt_inner_clusters_fingerprints_maxscore &a, const utils::dt_inner_clusters_fingerprints_maxscore &b) {
        return a.score > b.score;
    });
    return filteredPopCountBinsWithMaxScore;
}

std::tuple<std::vector<std::tuple<std::string, float>>, uint64_t> FPSearchEngine::search(const std::string& query, float threshold, int limits) {
    uint64_t *queryCFp = prepareQuery(query);
    std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore = filterPopcountBins(queryCFp[_fpStore->_CFPPopCountIndex], threshold);
    std::vector<std::tuple<std::string, float>> results;
    uint64_t num_sim_ops = 0;

    (this->*_normal_search)(filteredPopCountBinsWithMaxScore, queryCFp, threshold, limits, results, &num_sim_ops);

    //sort results by score
    std::sort(results.begin(), results.end(), [](const std::tuple<std::string, float>& a, const std::tuple<std::string, float>& b) {
        return std::get<1>(a) > std::get<1>(b);
    });

    // only return top results up to the specified limit
    if (results.size() > static_cast<size_t>(limits)) {
        results.resize(limits);
    }

    delete[] queryCFp;
    py::gil_scoped_acquire acquire;
    return std::make_tuple(results, num_sim_ops);
}

/*
void FPSearchEngine::_normal_search_memory_new(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& popCountBinsWithMaxScore, 
                                    uint64_t *queryCFp, float threshold, int limits,
                                    std::vector<std::tuple<std::string, float>> &results, uint64_t *num_sim_ops) {
    uint64_t commonPopCountThreshold = 0;
    float coeff, max_coeff = 0;
    uint64_t common_popcnt = 0;
    int num_hits = 0;
    constexpr uint64_t kUncomputedPopcnt = std::numeric_limits<uint64_t>::max();
    float running_threshold = 0;

    // cache for common popcount of query with cluster fps for each popcount bin to avoid 
    // repeated computation of common popcount for same cluster fps when iterating multiple
    // times over the same popcount bin for different thresholds (since clusters that do 
    // not qualify the threshold in one iteration might qualify in another iteration with 
    // lower threshold)
    std::vector<std::vector<uint64_t>> common_popcnt_clusters_cache;
    common_popcnt_clusters_cache.reserve(popCountBinsWithMaxScore.size());
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        common_popcnt_clusters_cache.emplace_back(
            popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters,
            kUncomputedPopcnt);
    }

    // build vector of thresholds to iterate over starting from 1.0 and decreasing by 0.1 down 
    // to the input threshold (inclusive)
    constexpr float kStep = 0.1f;
    std::vector<float> threshold_vec;
    threshold_vec.reserve(11);
    threshold_vec.push_back(1.0f);

    // Integer tenths avoid repeated float subtraction and reallocations.
    const int min_tenth = std::max(0, static_cast<int>(std::ceil(threshold * 10.0f)));
    for (int tenth = 9; tenth >= min_tenth; --tenth) {
        threshold_vec.push_back(static_cast<float>(tenth) * kStep);
    }

    // to keep track of clusters already processed for each popcount bin across different 
    // thresholds iterations
    std::vector<std::vector<uint8_t>> clusters_done;
    clusters_done.reserve(popCountBinsWithMaxScore.size());
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        clusters_done.emplace_back(
            popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters,
            static_cast<uint8_t>(0));
    }

    // iterating over thresholds starting from 1.0 and decreasing down to the input threshold (inclusive)
    for (size_t t=0; t < threshold_vec.size() && num_hits < limits; t++) {
        running_threshold = threshold_vec[t];

        //std::cout <<"running_threshold: " << running_threshold << ", num_hits: " << num_hits << std::endl;

        // iterating over popcount bins in order of their max score with query starting from highest until 
        // the max score is higher than the current threshold and required hits are not found yet
        for( size_t i=0; i < popCountBinsWithMaxScore.size() && running_threshold <= popCountBinsWithMaxScore[i].score && num_hits < limits; i++ ) {

            utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = popCountBinsWithMaxScore[i].inner_clusters_fingerprints;
            commonPopCountThreshold = (uint64_t) ceil(running_threshold * popCountBinsWithMaxScore[i].maxPopCountWithQuery);

            //std::cout<< "popCount: " << inner_clusters_fingerprints.popCount << ", maxPopCountWithQuery: " << popCountBinsWithMaxScore[i].maxPopCountWithQuery 
            //        << ", commonPopCountThreshold: " << commonPopCountThreshold 
            //        << ", maxScore: " << popCountBinsWithMaxScore[i].score << std::endl;

            uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
            uint64_t *fp_ptr = inner_clusters_fingerprints.fp;
            uint64_t inner_start = 0;
            std::vector<uint64_t>& bin_common_popcnt_cache = common_popcnt_clusters_cache[i];
            for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) {

                common_popcnt = bin_common_popcnt_cache[cid];
                if (common_popcnt == kUncomputedPopcnt) { // if common popcount for this cluster fps is not computed yet, compute and cache it
                    common_popcnt = bitwise_and_popcount(clusterFp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
                    bin_common_popcnt_cache[cid] = common_popcnt;
                    (*num_sim_ops)++;
                }
            
                // if common popcount of query with cluster fps is higher than the threshold and 
                // cluster is not processed yet in previous iterations for higher thresholds, then
                // process the cluster for potential hits, otherwise skip the cluster
                if (common_popcnt >= commonPopCountThreshold && !clusters_done[i][cid]) {
                    clusters_done[i][cid] = 1;
                    uint64_t inner_end = clusterFp_ptr[0];
                    for (auto fp_idx = inner_start; fp_idx < inner_end; fp_idx += _CFPSize, fp_ptr += _CFPSize) {

                        common_popcnt = bitwise_and_popcount(fp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
                        coeff =  TanimotoCoeff(common_popcnt, queryCFp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);
                        if (coeff >= threshold) {
                            results.push_back(std::make_tuple(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength), coeff));

                            // counting running hits for early stopping of bins loop
                            if (coeff >= popCountBinsWithMaxScore[i].score) num_hits++;
                        }
                        (*num_sim_ops)++;
                    }
                } else {
                    fp_ptr += clusterFp_ptr[0] - inner_start;
                }
                inner_start = clusterFp_ptr[0];
            }
        }

        //re-count hits
        num_hits = 0;
        for (auto &r : results) {
            if (std::get<1>(r) >= running_threshold) {
                num_hits++;
            }
        }
    }
}
*/

void FPSearchEngine::_normal_search_memory(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& popCountBinsWithMaxScore, 
                                    uint64_t *queryCFp, float threshold, int limits,
                                    std::vector<std::tuple<std::string, float>> &results, uint64_t *num_sim_ops) {
    uint64_t commonPopCountThreshold = 0;
    float coeff, max_coeff = 0;
    uint64_t common_popcnt = 0;
    
    for( auto inner_clusters_fingerprints_maxScore : popCountBinsWithMaxScore ) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if hits count is equal to the required limit
        if (max_coeff >= maxScore) {
            int hits = 0;
            for (auto &r : results) {
                if (std::get<1>(r) >= maxScore) {
                    hits++;
                }
            }
            if (hits >= limits) break;
        }

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queryCFp[_CFPPopCountIndex])); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t *fp_ptr = inner_clusters_fingerprints.fp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) {
            common_popcnt = bitwise_and_popcount(clusterFp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
            (*num_sim_ops)++;

            if (common_popcnt >= commonPopCountThreshold) {
                uint64_t inner_end = clusterFp_ptr[0];
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    common_popcnt = bitwise_and_popcount(fp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);;

                    if (common_popcnt >= commonPopCountThreshold) {
                        coeff =  TanimotoCoeff(common_popcnt, queryCFp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);
                        if (coeff >= threshold) {
                            results.push_back(std::make_tuple(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength), coeff));
                            if (coeff > max_coeff) max_coeff = coeff;
                        }
                    }
                    (*num_sim_ops)++;
                }
            } else {
                fp_ptr += clusterFp_ptr[0] - inner_start;
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}

void FPSearchEngine::_normal_search_disk(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& popCountBinsWithMaxScore, 
                                    uint64_t *queryCFp, float threshold, int limits,
                                    std::vector<std::tuple<std::string, float>> &results, uint64_t *num_sim_ops) {
    uint64_t commonPopCountThreshold = 0;
    float coeff, max_coeff = 0;
    uint64_t common_popcnt = 0;
    constexpr uint64_t kUncomputedPopcnt = std::numeric_limits<uint64_t>::max();

    std::vector<std::vector<uint64_t>> common_popcnt_clusters_cache;
    common_popcnt_clusters_cache.reserve(popCountBinsWithMaxScore.size());
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        common_popcnt_clusters_cache.emplace_back(
            popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters,
            kUncomputedPopcnt);
    }

    int num_hits = 0;
    //int binsLoopOver = 0;
    constexpr float kStep = 0.1f;
    std::vector<float> threshold_vec;
    threshold_vec.reserve(11);
    threshold_vec.push_back(1.0f);

    // Integer tenths avoid repeated float subtraction and reallocations.
    const int min_tenth = std::max(0, static_cast<int>(std::ceil(threshold * 10.0f)));
    for (int tenth = 9; tenth >= min_tenth; --tenth) {
        threshold_vec.push_back(static_cast<float>(tenth) * kStep);
    }

    std::vector<std::vector<uint8_t>> clusters_done;
    clusters_done.reserve(popCountBinsWithMaxScore.size());
    for (size_t i = 0; i < popCountBinsWithMaxScore.size(); ++i) {
        clusters_done.emplace_back(
            popCountBinsWithMaxScore[i].inner_clusters_fingerprints.num_clusters,
            static_cast<uint8_t>(0));
    }

    float running_threshold = 0;
    for (size_t t=0; t < threshold_vec.size() && num_hits < limits; t++) {
        running_threshold = threshold_vec[t];

        //std::cout <<"running_threshold: " << running_threshold << ", num_hits: " << num_hits << std::endl;

        for( size_t i=0; i < popCountBinsWithMaxScore.size() && running_threshold <= popCountBinsWithMaxScore[i].score && num_hits < limits; i++ ) {

            utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = popCountBinsWithMaxScore[i].inner_clusters_fingerprints;
            commonPopCountThreshold = (uint64_t) ceil(running_threshold * popCountBinsWithMaxScore[i].maxPopCountWithQuery);

            //std::cout<< "popCount: " << inner_clusters_fingerprints.popCount << ", maxPopCountWithQuery: " << popCountBinsWithMaxScore[i].maxPopCountWithQuery 
            //        << ", commonPopCountThreshold: " << commonPopCountThreshold 
            //        << ", maxScore: " << popCountBinsWithMaxScore[i].score << std::endl;

            uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
            uint64_t inner_start = 0;
            std::vector<uint64_t>& bin_common_popcnt_cache = common_popcnt_clusters_cache[i];
            for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) {

                common_popcnt = bin_common_popcnt_cache[cid];
                if (common_popcnt == kUncomputedPopcnt) {
                    common_popcnt = bitwise_and_popcount(clusterFp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
                    bin_common_popcnt_cache[cid] = common_popcnt;
                    (*num_sim_ops)++;
                }
            
                if (common_popcnt >= commonPopCountThreshold && !clusters_done[i][cid]) {
                    uint64_t inner_end = clusterFp_ptr[0];
                    uint64_t *fp_ptr = _fpStore->getFPsForCluster(inner_clusters_fingerprints.popCount, inner_start, inner_end); // read fps for this cluster from disk
                    for (auto fp_idx = inner_start; fp_idx < inner_end; fp_idx += _CFPSize, fp_ptr += _CFPSize) {

                        common_popcnt = bitwise_and_popcount(fp_ptr+_molIdOffset, queryCFp+_molIdOffset, _fpSize);
                        coeff =  TanimotoCoeff(common_popcnt, queryCFp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);
                        if (coeff >= threshold) {
                            results.push_back(std::make_tuple(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength), coeff));

                            // counting running hits for early stopping of bins loop
                            if (coeff >= popCountBinsWithMaxScore[i].score) num_hits++;
                        }
                        (*num_sim_ops)++;
                    }
                    free(fp_ptr - (inner_end - inner_start)); // free memory allocated for fps read from disk
                }
                inner_start = clusterFp_ptr[0];
            }
        }

        //re-count hits
        num_hits = 0;
        for (auto &r : results) {
            if (std::get<1>(r) >= running_threshold) {
                num_hits++;
            }
        }
    }
}

std::vector<std::vector<std::tuple<std::string, float>>> FPSearchEngine::batchSearch(const std::vector<std::string>& queries, float threshold, int limits) {
    // build query Compact Fingerprints
    std::vector<uint64_t*> queriesCFp;
    queriesCFp.reserve(queries.size());
    for (size_t i=0; i < queries.size(); i++) {
        uint64_t *queryCFp = prepareQuery(queries[i]);
        queryCFp[0] = i;
        queriesCFp.push_back(queryCFp);
    }

    //sort queries by popcount
    std::sort(queriesCFp.begin(), queriesCFp.end(), [this](uint64_t *a, uint64_t *b) {
        return a[_CFPPopCountIndex] > b[_CFPPopCountIndex];
    });

    //prepare batch data by sorting queries by popcount
    uint64_t previous_popcount = 0;
    std::vector<utils::dt_batch_data> batch_data;
    for (size_t i=0; i < queriesCFp.size(); i++) {
        if (queriesCFp[i][_CFPPopCountIndex] != previous_popcount) { // create new batch data
            utils::dt_batch_data new_batch_data;         
            std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filteredPopCountBinsWithMaxScore = filterPopcountBins(queriesCFp[i][_CFPPopCountIndex], threshold);
            new_batch_data.filteredPopCountBinsWithMaxScore = filteredPopCountBinsWithMaxScore;
            new_batch_data.popCount = queriesCFp[i][_CFPPopCountIndex];
            batch_data.push_back(new_batch_data);
        }
        // create new query data
        utils::dt_batch_query_data new_query_data;
        new_query_data.cfp = queriesCFp[i];

        // add new query data to batch data
        utils::dt_batch_query_data *qdata = batch_data.at(batch_data.size()-1).qdata;
        int qsize = batch_data.at(batch_data.size()-1).qsize;
        qdata = (utils::dt_batch_query_data*) realloc(qdata, sizeof(utils::dt_batch_query_data) * (qsize + 1)); // memory allocation
        qdata[qsize] = new_query_data;
        batch_data.at(batch_data.size()-1).qdata = qdata;
        batch_data.at(batch_data.size()-1).qsize = qsize + 1;

        previous_popcount = queriesCFp[i][_CFPPopCountIndex];
    }

    std::vector<std::vector<std::tuple<std::string, float>>> finalResults;
    finalResults.resize(queries.size());
    for(auto bdata : batch_data) {
        // perform batch search
        (this->*_batch_search)(bdata, threshold, limits);

        // process and assign results
        for(auto q=0; q < bdata.qsize; q++) {
            // prepare results
            std::vector<std::tuple<std::string, float>> results;
            results.reserve(bdata.qdata[q].results_size);
            for(auto r=0; r < bdata.qdata[q].results_size; r++) {
                results.push_back(std::make_tuple(*bdata.qdata[q].results[r].id, bdata.qdata[q].results[r].score));
            }
            
            // sort results by score
            std::sort(results.begin(), results.end(), [](const std::tuple<std::string, float>& a, const std::tuple<std::string, float>& b) {
                return std::get<1>(a) > std::get<1>(b);
            });

            // assign results to finalResults
            finalResults[bdata.qdata[q].cfp[0]] = results;
        }
    }

    //release memory for batch_data
    for(auto bdata : batch_data) {
        for(auto q=0; q < bdata.qsize; q++) {
            for(auto r=0; r < bdata.qdata[q].results_size; r++)
                delete bdata.qdata[q].results[r].id; // id created through new operator
            free(bdata.qdata[q].results); // by realloc
        }
        free(bdata.qdata); // by realloc
    }

    // remove memory for queryCFp
    for (size_t i=0; i < queriesCFp.size(); i++) {
        delete[] queriesCFp[i];
    }

    py::gil_scoped_acquire acquire;
    return finalResults;
}

void FPSearchEngine::_batch_search_memory(utils::dt_batch_data &batch_data, float threshold, int limits) {
    uint64_t commonPopCountThreshold = 0;
    uint64_t queriesPopcount = batch_data.popCount;
    float coeff = 0;
    uint64_t max_common_popcnt = 0;
    uint64_t common_popcnt = 0;
    bool all_done = false;
    utils::dt_batch_query_data *query_data = batch_data.qdata;
    
    for( auto &inner_clusters_fingerprints_maxScore : batch_data.filteredPopCountBinsWithMaxScore) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if all queries are done within requested limits
        all_done = true;
        query_data = batch_data.qdata;
        for(auto q=0; q < batch_data.qsize; q++, query_data++) {
            if (!query_data->done && query_data->max_coeff >= maxScore) {
                int hits = 0;
                for(auto r=0; r < query_data->results_size; r++) {
                    if (query_data->results[r].score >= maxScore) hits++;
                }
                if (hits >= limits) query_data->done = true;
            }
            all_done &= query_data->done;
        }
        if (all_done) break;

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queriesPopcount)); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t *fp_ptr = inner_clusters_fingerprints.fp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) { //searching through clusters

            // check if any query need to be search in current cluster
            max_common_popcnt = 0;
            query_data = batch_data.qdata;
            for(auto q=0; q < batch_data.qsize; q++, query_data++) {
                if (!query_data->done) {
                    common_popcnt = bitwise_and_popcount(clusterFp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                    //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                    //    common_popcnt += popcntll(clusterFp_ptr[j] & query_data->cfp[j]);
                    //}
                    if (common_popcnt > max_common_popcnt) max_common_popcnt = common_popcnt;
                }
            }

            // queries are searched in current cluster
            if (max_common_popcnt >= commonPopCountThreshold) { // potentially hit could be found in current cluster for at least one query
                uint64_t inner_end = clusterFp_ptr[0];
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    query_data = batch_data.qdata;
                    for(auto q=0; q < batch_data.qsize; q++, query_data++) { // similarity of each query with each fingerprint of current cluster
                        if (!query_data->done) {
                            common_popcnt = bitwise_and_popcount(fp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                            //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                            //    common_popcnt += popcntll(fp_ptr[j] & query_data->cfp[j]);
                            //}

                            if (common_popcnt >= commonPopCountThreshold) { // potential hit found
                                coeff =  TanimotoCoeff(common_popcnt, query_data->cfp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);

                                if (coeff >= threshold) { // exact hit found, add to results
                                    utils::dt_result *results = query_data->results;
                                    results = (utils::dt_result*) realloc(results, sizeof(utils::dt_result) * (query_data->results_size + 1));
                                    std::string *id = new std::string(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength));
                                    results[query_data->results_size] = {id, coeff};
                                    query_data->results = results;
                                    query_data->results_size += 1;
                                    if (coeff > query_data->max_coeff) query_data->max_coeff = coeff;
                                }
                            }
                        }
                    }
                }
            } else {
                fp_ptr += clusterFp_ptr[0] - inner_start;
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}

void FPSearchEngine::_batch_search_disk(utils::dt_batch_data &batch_data, float threshold, int limits) {
    uint64_t commonPopCountThreshold = 0;
    uint64_t queriesPopcount = batch_data.popCount;
    float coeff = 0;
    uint64_t max_common_popcnt = 0;
    uint64_t common_popcnt = 0;
    bool all_done = false;
    utils::dt_batch_query_data *query_data = batch_data.qdata;
    
    for( auto &inner_clusters_fingerprints_maxScore : batch_data.filteredPopCountBinsWithMaxScore) {
        float maxScore = inner_clusters_fingerprints_maxScore.score;

        // check if all queries are done within requested limits
        all_done = true;
        query_data = batch_data.qdata;
        for(auto q=0; q < batch_data.qsize; q++, query_data++) {
            if (!query_data->done && query_data->max_coeff >= maxScore) {
                int hits = 0;
                for(auto r=0; r < query_data->results_size; r++) {
                    if (query_data->results[r].score >= maxScore) hits++;
                }
                if (hits >= limits) query_data->done = true;
            }
            all_done &= query_data->done;
        }
        if (all_done) break;

        utils::dt_inner_clusters_fingerprints inner_clusters_fingerprints = inner_clusters_fingerprints_maxScore.inner_clusters_fingerprints;
        commonPopCountThreshold = (uint64_t) ceil(threshold * std::max(inner_clusters_fingerprints.popCount, (int)queriesPopcount)); // get the common popcount threshold for this bin

        uint64_t *clusterFp_ptr = inner_clusters_fingerprints.clusterFp;
        uint64_t inner_start = 0;
        for(size_t cid=0; cid < inner_clusters_fingerprints.num_clusters; cid++, clusterFp_ptr += _CFPSize) { //searching through clusters

            // check if any query need to be search in current cluster
            max_common_popcnt = 0;
            query_data = batch_data.qdata;
            for(auto q=0; q < batch_data.qsize; q++, query_data++) {
                if (!query_data->done) {
                    common_popcnt = bitwise_and_popcount(clusterFp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                    //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                    //    common_popcnt += popcntll(clusterFp_ptr[j] & query_data->cfp[j]);
                    //}
                    if (common_popcnt > max_common_popcnt) max_common_popcnt = common_popcnt;
                }
            }

            // queries are searched in current cluster
            if (max_common_popcnt >= commonPopCountThreshold) { // potentially hit could be found in current cluster for at least one query
                uint64_t inner_end = clusterFp_ptr[0];
                uint64_t *fp_ptr = _fpStore->getFPsForCluster(inner_clusters_fingerprints.popCount, inner_start, inner_end); // read fps for this cluster from disk
                for (auto i = inner_start; i < inner_end; i+=_CFPSize, fp_ptr += _CFPSize) {

                    query_data = batch_data.qdata;
                    for(auto q=0; q < batch_data.qsize; q++, query_data++) { // similarity of each query with each fingerprint of current cluster
                        if (!query_data->done) {
                            common_popcnt = bitwise_and_popcount(fp_ptr + _molIdOffset, query_data->cfp + _molIdOffset, _fpSize);
                            //for (auto j = _molIdOffset; j < _fpEndIndex; j++) {
                            //    common_popcnt += popcntll(fp_ptr[j] & query_data->cfp[j]);
                            //}

                            if (common_popcnt >= commonPopCountThreshold) { // potential hit found
                                coeff =  TanimotoCoeff(common_popcnt, query_data->cfp[_CFPPopCountIndex], fp_ptr[_CFPPopCountIndex], _div_lookup_table);

                                if (coeff >= threshold) { // exact hit found, add to results
                                    utils::dt_result *results = query_data->results;
                                    results = (utils::dt_result*) realloc(results, sizeof(utils::dt_result) * (query_data->results_size + 1));
                                    std::string *id = new std::string(utils::getMolIdFromCompactFPArray(fp_ptr, _molIdMaxLength));
                                    results[query_data->results_size] = {id, coeff};
                                    query_data->results = results;
                                    query_data->results_size += 1;
                                    if (coeff > query_data->max_coeff) query_data->max_coeff = coeff;
                                }
                            }
                        }
                    }
                }
                free(fp_ptr - (inner_end - inner_start)); // free memory allocated for fps read from disk
            }
            inner_start = clusterFp_ptr[0];
        }
    }
}
