// lib/Utils.groovy

class Utils {

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