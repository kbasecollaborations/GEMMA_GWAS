/*
A KBase module: GEMMA_GWAS
*/

module GEMMA_GWAS {
    /*
        KBase style object reference X/Y/Z to a KBaseMatrices.TraitMatrix structure
            @id ws KBaseMatrices.TraitMatrix
    */
    typedef string trait_ref;

    /*
        KBase style object reference X/Y/Z to a
            @id ws KBaseGwasData.Variations
    */
    typedef string var_ref;

	typedef structure {
	    string workspace_name;
		trait_ref trait_matrix;
		var_ref variation;
	} GemmaGwasInput;

    typedef structure {
        string report_name;
        string report_ref;
    } GwasResults;

    funcdef run_gemma_association(GemmaGwasInput params) returns (GwasResults output) authentication required;
};
