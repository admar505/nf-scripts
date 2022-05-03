Tim Diels @timdiels Nov 06 2017 06:30
// Parse species.txt
Channel
    .fromPath("$params.input/species.txt")
    .splitCsv(sep: "\t", header: true)
    .map {
        def specie = [
            name: it.name,
            id: it.id,
            proteome: Paths.get(it.proteome),
            longName: it.long_name,
            orthofinder: it.orthofinder.toBoolean(),
            phyml: it.phyml.toBoolean(),
            interproscan: it.interproscan.toBoolean(),
            blast2go: it.blast2go.toBoolean(),
            clime: it.clime.toBoolean(),
        ]

        // Assert proteome is an absolute path
        if (!specie.proteome.isAbsolute()) {
            error(
                "$specie.name's proteome must be an absolute path, got: $specie.proteome"
            )
        }

        // Assert OrthoFinder is set when PhyML is set
        if (!specie.orthofinder && specie.phyml) {
            error(
                "$specie.name has PhyML set, but not OrthoFinder. PhyML requires " +
                "OrthoFinder. Please set OrthoFinder as well or unset PhyML."
            )
        }

        // Assert PhyML is set when CLIME is set
        if (!specie.phyml && specie.clime) {
            error(
                "$specie.name has CLIME set, but not PhyML. CLIME requires " +
                "PhyML. Please set PhyML as well or unset CLIME."
            )
        }

        //
        return specie
    }
    .set { species }
