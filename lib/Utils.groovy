// lib/Utils.groovy

import java.nio.file.Path
import nextflow.Global
import nextflow.Nextflow

class Utils {

    // Copy an input file (URL, cloud, or local path) into the work dir once and reuse that copy.
    // Tasks then stage it from the work dir instead of the original source each time, which
    // avoids rate limits (e.g. GitHub) on large runs. On -resume the existing copy is reused.
    static Path stageOnce(String src) {
        def dest = Global.session.workDir.resolve("assets/${src.md5()}/${Nextflow.file(src).name}")
        if( !dest.exists() ) Nextflow.file(src, checkIfExists: true).copyTo(dest)
        return dest
    }

    static String getRefName(filename, String sampleName) {
        // Convert to string and strip any path
        def name = filename.getName()

        // Strip sampleName prefix: "<sampleName>_"
        if (sampleName && name.startsWith(sampleName + "_")) {
            name = name.substring((sampleName + "_").length())
        }

        // Handle .gz first
        if (name.endsWith(".gz")) {
            name = name[0..-4]
        }

        // Now strip FASTA-style extensions
        def fastaExts = [".fa", ".fna", ".fasta"]
        def ext = fastaExts.find { name.toLowerCase().endsWith(it) }
        if (ext) {
            name = name[0..-(ext.size() + 1)]
        }

        return name
    }
}