
process guide_assignment_cleanser {

    input:
        path mudata_input
        val threshold
    output:
        path "cleanser_mudata_output.h5mu", emit: guide_assignment_mudata_output

    script:
        def thresh_opt = threshold ? "-t ${threshold}" : ""
        """
        export CMDSTAN=/root/.cmdstan/cmdstan-2.36.0
        export PATH=\$PATH:\$CMDSTAN/bin
        
        cleanser -i ${mudata_input} --posteriors-output cleanser_mudata_output.h5mu --modality guide --capture-method capture_method --output-layer guide_assignment ${thresh_opt}
        """
}