"""
Main search-engine for UFFPSim.

Author: Rajendra Kumar

"""

from typing import List, Tuple
import json

from . import uffpsimLib
from .fingerprints import FPCalculator, load_molecule

class UFFPSimSearchEngine(uffpsimLib.FPSearchEngineBase):
    """A class representing a fingerprint-based search engine for UFFPSim.

    Attributes
    ----------
    db_file : str
        The file path of the database containing the molecular fingerprints and their corresponding identifiers.
    running : int
        A counter to keep track of the number of search operations currently in progress.
    fp_calculator : FPCalculator
        An instance of the FPCalculator class for calculating molecular fingerprints.

    Example
    -------
    >>> search_engine = UFFPSimSearchEngine("database.db")
    >>> result = search_engine.search("Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl", 0.7, 1)
    >>> print(result)
    >>> # Batch search operation:
    >>> query_smiles = ["Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl", "Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1"]
    >>> batch_result = search_engine.batch_search(query_smiles, 0.7, 1)
    >>> [print(q, r) for q, r in zip(query_smiles, batch_result)]
    """
    def __init__(self, db_file: str, mode: str = "memory"):
        """
        Initializes the UFFPSimSearchEngine with the given database file and initializes the fingerprint calculator.

        Parameters
        ----------
        db_file : str
            The file path of the database containing the molecular fingerprints and their corresponding identifiers.
        mode : str, optional
            The mode of operation for the search engine. Can be "memory" or "disk". Default is "memory".

        """
        super().__init__(db_file, mode)
        self.running = 0
        fp_input_arguments = json.loads(self.fp_store.fp_params_json)
        self.fp_calculator = FPCalculator(fp_input_arguments["fp_type"], fp_input_arguments["fp_params"])

    def search(self, mol_data: str, threshold: float = 0.2, limit_by: int = 10) -> Tuple[List[Tuple[str, float]],int]:
        """
        Performs a single search operation using the given molecular data and returns a list of tuples containing the identifiers and similarity scores of the top-k matches.

        Parameters
        ----------
        mol_data : str
            The molecular data for which the search operation is performed.
        threshold : float, optional
            The similarity threshold for filtering the search results. Default is 0.2.
        limit_by : int, optional
            The maximum number of search results to return. Default is 10.

        Returns
        -------
        None | List[Tuple[str, float]]
            A list of tuples containing the identifiers and similarity scores of the top-k matches.
            If an error occurs while loading a molecule, it is logged and the corresponding result is set to None.

        """
        self.running += 1
        mol, smiles = load_molecule(mol_data)
        if mol is not None:
            fp = self.fp_calculator(mol)
            self.running -= 1
            return self._search(fp, threshold, limit_by)
        else:
            print(f"Error loading molecule: {mol_data}")
            return None
    
    def batch_search(self, mols_data: [str], threshold: float = 0.2, limit_by: int = 10) -> List[Tuple[str, float]]:
        """
        Performs a batch search operation using the given list of molecular data and returns a list of tuples containing the identifiers and similarity scores of the top-k matches for each molecule.

        Parameters
        ----------
        mols_data : List[str]
            The list of molecular data for which the batch search operation is performed.
        threshold : float, optional
            The similarity threshold for filtering the search results. Default is 0.2.
        limit_by : int, optional
            The maximum number of search results to return for each molecule. Default is 10.

        Returns
        -------
        List[None | List[Tuple[str, float]]]
            A list of tuples containing the identifiers and similarity scores of the top-k matches for each molecule.
            If an error occurs while loading a molecule, it is logged and the corresponding result is set to None.

        """
        self.running += 1
        fingerprints = []
        searched_indices = []
        for idx, mol_data in enumerate(mols_data):
            mol, smiles = load_molecule(mol_data)
            if mol is not None:
                fp = self.fp_calculator(mol)
                fingerprints.append(fp)
                searched_indices.append(idx)
            else:
                print(f"Error loading molecule: {mol_data}")

        raw_result = self._batch_search(fingerprints, threshold, limit_by)
        final_results = [None] * len(mols_data)
        for idx, result in zip(searched_indices, raw_result):
            final_results[idx] = result
        self.running -= 1

        return final_results

    def getCompactFingerPrintArray(self, mol_data):
        """
        Calculates and returns the compact fingerprint array for the given molecular data.

        Parameters
        ----------
        mol_data : str
            The molecular data for which the compact fingerprint array is calculated.

        Returns
        -------
        List[int]
            The compact fingerprint array for the given molecular data.

        """
        fp = self.fp_calculator(load_molecule(mol_data)[0])
        return uffpsimLib.getCompactFingerPrintArray(fp)
    
    def get_smiles_for_id(self, id):
        """
        Retrieves the SMILES representation of a molecule using its identifier.

        Parameters
        ----------
        id : str
            The identifier of the molecule.

        Returns
        -------
        str
            The SMILES representation of the molecule.

        """
        return self.fp_store.get_smiles_for_id(id)
    
    def build_mol_id_to_index_map(self):
        """
        Builds a mapping of molecule identifiers to their corresponding indices in the database.

        Returns
        -------
        None

        """
        return self.fp_store.build_mol_id_to_index_map()