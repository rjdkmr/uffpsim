/**
 * @brief This class represents a search engine for molecular fingerprints.
 * @author Rajendra Kumar
 */

#include <iostream>
#include <cstdint>
#include <map>
#include <vector>
#include <tuple>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

#include "fpstore.h"

#define TanimotoCoeff(common, query, target, lookup_table) common*lookup_table[query+target-common]

class FPSearchEngine {
    private:

        using NormalSearchFn = void (FPSearchEngine::*)(
            const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>&,
            uint64_t *,
            float,
            int,
            std::vector<std::tuple<std::string, float>>&,
            uint64_t *);

        NormalSearchFn _normal_search = nullptr;

        using BatchSearchFn = void (FPSearchEngine::*)(
            utils::dt_batch_data &,
            float,
            int);

        BatchSearchFn _batch_search = nullptr;

        /**
         * The _normal_search method performs a single-query search on the fingerprint store using the Tanimoto coefficient.
         * It compares the query fingerprint with the fingerprints of the inner clusters and populates the results vector with 
         * the molecule IDs and similarity scores.
         *
         * @param inner_clusters_fingerprints_maxscore: A vector of inner clusters' fingerprints popcount-wise and their maximum scores.
         * @param queryCFp: A pointer to the query fingerprint.
         * @param threshold: The similarity threshold for filtering the results.
         * @param limits: The maximum number of results to return.
         * @param results: A reference to a vector of tuples, where each tuple contains the molecule ID and similarity score.
         *
         * Note: The _normal_search method assumes that the inner_clusters_fingerprints_maxscore vector is sorted by the maximum score in descending order.
         */
        void _normal_search_memory(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& inner_clusters_fingerprints_maxscore, uint64_t *queryCFp, 
                            float threshold, int limits, std::vector<std::tuple<std::string, float>> &results, uint64_t *num_sim_ops);

        void _normal_search_disk(const std::vector<utils::dt_inner_clusters_fingerprints_maxscore>& inner_clusters_fingerprints_maxscore, uint64_t *queryCFp, 
                    float threshold, int limits, std::vector<std::tuple<std::string, float>> &results, uint64_t *num_sim_ops);

        /**
         * The prepareQuery method prepares the query fingerprint for searching.
         * It converts the query fingerprint string into a 64-bit unsigned integer representation.
         *
         * @param query: The query fingerprint string.
         *
         * @return: A pointer to the prepared query fingerprint.
         *
         * Note: The prepareQuery method assumes that the query fingerprint string is in a valid format and 
         *       can be converted into a 64-bit unsigned integer.
         */
        uint64_t* prepareQuery(const std::string& query);

        /**
         * The filterPopcountBins method filters and sort the the popcount-bins using the similarity threshold.
         * It select popcount with their inner-clusters falls within a certain range and whose similarity score with 
         * the query fingerprint is above the threshold.
         *
         * @param queryPopcount: The popcount of the query fingerprint.
         * @param threshold: The similarity threshold for filtering the results.
         *
         * @return: A vector of inner clusters' fingerprints popcount-wise and their maximum scores that satisfy the filtering criteria.
         */
        std::vector<utils::dt_inner_clusters_fingerprints_maxscore> filterPopcountBins(uint64_t queryPopcount, float threshold);
        
        /**
         * The _batch_search method performs a batch search on the fingerprint store using the Tanimoto coefficient.
         * It compares multiple query fingerprints with the fingerprints of the inner clusters and populates the results vector with 
         * the molecule IDs and similarity scores for each query.
         *
         * @param batch_data: A reference to a structure containing the batch data, including query fingerprints and their corresponding molecule IDs.
         * @param threshold: The similarity threshold for filtering the results.
         * @param limits: The maximum number of results to return for each query.
         *
         * Note: The _batch_search method assumes that the batch_data structure is properly initialized with valid query fingerprints and molecule IDs.
         *       It also assumes that the inner_clusters_fingerprints_maxscore vector is sorted by the maximum score in descending order.
         */
        void _batch_search_memory(utils::dt_batch_data &batch_data, float threshold, int limits);

        void _batch_search_disk(utils::dt_batch_data &batch_data, float threshold, int limits);

    public:

        /**
         * The FPSearchEngine class is responsible for performing molecular similarity searches using fingerprints.
         * It utilizes a pre-built fingerprint store and provides methods for single-query and batch searches.
         *
         * @param filename: The path to the pre-built fingerprint store file.
         *
         * Note: The FPSearchEngine class assumes that the fingerprint store file is properly formatted and contains valid fingerprint data.
         */
        FPSearchEngine(const std::string& filename, std::string mode = "memory");
        ~FPSearchEngine();

        std::string _mode = "memory"; // default search mode is memory, it can be set to "disk" for disk-based search
        int _molIdMaxLength;
        int _fpSize = 0;
        int _CFPSize;
        int _molIdOffset;
        int _CFPPopCountIndex;
        int _fpEndIndex;
        std::vector<float> _div_lookup_table;

        FingerprintStore *_fpStore; // pointer to the fingerprint store

        /**
         * The close method closes the fingerprint store file and releases any associated resources.
         *
         * Note: The close method should be called when the FPSearchEngine object is no longer needed to free up system resources.
         */
        void close();

        /**
         * The getNumberOfMolecules method returns the total number of molecules in the fingerprint store.
         *
         * @return: The number of molecules in the fingerprint store.
         *
         */
        uint64_t getNumberOfMolecules() {return _fpStore->getNumberOfMolecules(); };
        
        /**
         * The search method performs a single-query molecular similarity search using the Tanimoto coefficient.
         *
         * @param fp_string: The query fingerprint string.
         * @param threshold: The similarity threshold for filtering the results.
         * @param limits: The maximum number of results to return.
         *
         * @return: A vector of tuples, where each tuple contains the molecule ID and similarity score.
         *
         */
        std::tuple<std::vector<std::tuple<std::string, float>>, uint64_t> search(const std::string &fp_string, float threshold, int limits);
        
        /**
         * The batchSearch method performs a batch molecular similarity search using the Tanimoto coefficient.
         *
         * @param queries: A vector of query fingerprint strings.
         * @param threshold: The similarity threshold for filtering the results.
         * @param limits: The maximum number of results to return for each query.
         *
         * @return: A vector of vectors, where each inner vector contains tuples, and each tuple contains the molecule ID and similarity score for a query.
         *
         */
        std::vector<std::vector<std::tuple<std::string, float>>> batchSearch(const std::vector<std::string>& queries, float threshold, int limits);
};