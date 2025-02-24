
process CreateMuData {

  cache 'lenient'

  input:

  path adata_rna
  path adata_guide
  path guide_metadata
  path gtf_file
  val moi

  output:
  path "mudata.h5mu" , emit: mudata
  path "guide_concatenated_adata.h5ad", emit: adata_guide
  
  script:
        """
        pip install git+https://github.com/LucasSilvaFerreira/GTFProcessing.git
        create_mdata.py ${adata_rna} ${adata_guide} ${guide_metadata} ${gtf_file} ${moi}
        mv concatenated_adata.h5ad guide_concatenated_adata.h5ad
        """

}
