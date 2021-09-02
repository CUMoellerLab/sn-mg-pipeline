import numpy as np
import pandas as pd
import gzip
from os import path
from skbio import DistanceMatrix
from yaml import dump

def _validate_parameters(dm, num_prototypes, seedset=None):
    '''Validate the paramters for each algorithm.
    Parameters
    ----------
    dm: skbio.stats.distance.DistanceMatrix
        Pairwise distances for all elements in the full set S.
    num_prototypes: int
        Number of prototypes to select for distance matrix.
        Must be >= 2, since a single prototype is useless.
        Must be smaller than the number of elements in the distance matrix,
        otherwise no reduction is necessary.
    seedset: iterable of str
        A set of element IDs that are pre-selected as prototypes. Remaining
        prototypes are then recruited with the prototype selection algorithm.
        Warning: It will most likely violate the global objective function.
    Raises
    ------
    ValueError
        The number of prototypes to be found should be at least 2 and at most
        one element smaller than elements in the distance matrix. Otherwise, a
        ValueError is raised.
        The IDs in the seed set must be unique, and must be present in the
        distance matrix. Otherwise, a ValueError is raised.
        The size of the seed set must be smaller than the number of prototypes
        to be found. Otherwise, a ValueError is raised.
    '''
    if num_prototypes < 2:
        raise ValueError("'num_prototypes' must be >= 2, since a single "
                         "prototype is useless.")
    if num_prototypes >= dm.shape[0]:
        raise ValueError("'num_prototypes' must be smaller than the number of "
                         "elements in the distance matrix, otherwise no "
                         "reduction is necessary.")
    if seedset is not None:
        seeds = set(seedset)
        if len(seeds) < len(seedset):
            raise ValueError("There are duplicated IDs in 'seedset'.")
        if not seeds < set(dm.ids):  # test if set A is a subset of set B
            raise ValueError("'seedset' is not a subset of the element IDs in "
                             "the distance matrix.")
        if len(seeds) >= num_prototypes:
            raise ValueError("Size of 'seedset' must be smaller than the "
                             "number of prototypes to select.")

def prototype_selection_destructive_maxdist(dm, num_prototypes, seedset=None):
    '''Heuristically select k prototypes for given distance matrix.

       Prototype selection is NP-hard. This is an implementation of a greedy
       correctness heuristic: Start with the complete set and iteratively
       remove elements until the number of required prototypes is left.
       The decision which element shall be removed is based on the minimal
       distance sum this element has to all other.

    Parameters
    ----------
    dm: skbio.stats.distance.DistanceMatrix
        Pairwise distances for all elements in the full set S.
    num_prototypes: int
        Number of prototypes to select for distance matrix.
        Must be >= 2, since a single prototype is useless.
        Must be smaller than the number of elements in the distance matrix,
        otherwise no reduction is necessary.
    seedset: iterable of str
        A set of element IDs that are pre-selected as prototypes. Remaining
        prototypes are then recruited with the prototype selection algorithm.
        Warning: It will most likely violate the global objective function.

    Returns
    -------
    list of str
        A sequence holding selected prototypes, i.e. a sub-set of the
        IDs of the elements in the distance matrix.

    Raises
    ------
    ValueError
        The number of prototypes to be found should be at least 2 and at most
        one element smaller than elements in the distance matrix. Otherwise, a
        ValueError is raised.

    Notes
    -----
    Timing: %timeit -n 100 prototype_selection_destructive_maxdist(dm, 100)
            100 loops, best of 3: 2.1 s per loop
            where the dm holds 27,398 elements
    function signature with type annotation for future use with python >= 3.5:
    def prototype_selection_constructive_maxdist(dm: DistanceMatrix,
    num_prototypes: int, seedset: List[str]) -> List[str]:

    ...
    LICENSE
    From https://github.com/biocore/wol/tree/master/code/prototypeSelection

    Copyright (c) 2017--, WoL development team. All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the name of the copyright holder nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    '''
    _validate_parameters(dm, num_prototypes, seedset)

    # clever bookkeeping allows for significant speed-ups!

    # track the number of available elements
    numRemain = len(dm.ids)

    # distances from each element to all others
    currDists = dm.data.sum(axis=1)

    # a dirty hack to ensure that all elements of the seedset will be selected
    # last and thus make it into the resulting set
    maxVal = currDists.max()
    if seedset is not None:
        for e in seedset:
            currDists[dm.index(e)] = maxVal*2

    # the element to remove first is the one that has smallest distance to all
    # other. "Removing" works by tagging its distance-sum as infinity. Plus, we
    # decrease the number of available elements by one.
    minElmIdx = currDists.argmin()
    currDists[minElmIdx], numRemain = np.infty, numRemain-1

    # continue until only num_prototype elements are left
    while (numRemain > num_prototypes):
        # substract the distance to the removed element for all remaining
        # elements
        currDists -= dm.data[minElmIdx]
        # find the next element to be removed, again as the one that is
        # closest to all others
        minElmIdx = currDists.argmin()
        currDists[minElmIdx], numRemain = np.infty, numRemain-1

    # return a list of IDs of the surviving elements, which are the found
    # prototypes.
    return [dm.ids[idx]
            for idx, dist in enumerate(currDists)
            if dist != np.infty]


rule sourmash_sketch_reads:
    input:
        R1 = "output/filtered/nonhost/{sample}.1.fastq.gz",
        R2 = "output/filtered/nonhost/{sample}.2.fastq.gz"
    output:
        "output/sourmash/sketches/{sample}.sig"
    log:
        "output/logs/sourmash/sourmash_sketch_reads.{sample}.log"
    threads: 1
    conda: "../env/sourmash.yaml"
    params:
        k = config['params']['sourmash']['k'],
        scaled = config['params']['sourmash']['scaled'],
        extra = config['params']['sourmash']['extra']
    shell:
        """
        sourmash sketch dna \
        -p k={params.k},scaled={params.scaled}  \
        {params.extra} \
        -o {output} \
        --merge \
        {input} 2> {log} 1>&2
        """

rule sourmash_dm:
    input:
        expand(rules.sourmash_sketch_reads.output,
               sample=samples)
    output:
        dm = "output/sourmash/sourmash.dm",
        csv = "output/sourmash/sourmash.csv",
        labels = "output/sourmash/sourmash.dm.labels.txt"
    log:
        "output/logs/sourmash/sourmash_dm.log"
    threads: 1
    conda: "../env/sourmash.yaml"
    shell:
        """
        sourmash compare \
        --output {output.dm} \
        --csv {output.csv} \
        {input} 2> {log} 1>&2
        """

rule sourmash_plot:
    input:
        "output/sourmash/sourmash.dm"
    output:
        directory("output/sourmash/plots")
    log:
        "output/logs/sourmash/sourmash_plot.log"
    threads: 1
    conda: "../env/sourmash.yaml"
    shell:
        """
        sourmash plot --pdf --labels \
        --output-dir {output} \
        {input} 2> {log} 1>&2
        """

rule prototype_selection:
    input:
        dm = rules.sourmash_dm.output.csv,
        labels = rules.sourmash_dm.output.labels
    output:
        file = "output/sourmash/selected_prototypes.yaml"
    params:
        min_seqs = config['params']['prototypes']['min_seqs'],
        max_seqs = config['params']['prototypes']['max_seqs']
    log:
        "output/logs/sourmash/prototype_selection.log"
    threads: 1
    run:
        df = pd.read_csv(input[0], header=0, encoding= 'unicode_escape')
        df.index = df.columns

        # test file sizes

        pf_seqs = []
        for fp in df.columns:
            print(fp)
            with gzip.open(fp, 'rb') as f:
                for i, l in enumerate(f):
                    pass
            seqs = (i + 1) / 4
            print(seqs)
            if params['min_seqs'] <= seqs <= params['max_seqs']:
                pf_seqs.append(fp)

        df_filt = df.loc[pf_seqs, pf_seqs]

        labels = [os.path.basename(x) for x in pf_seqs]

        dm = DistanceMatrix(1 - df_filt.values)

        print("The imported distance matrix has "
              "{} elements.".format(len(labels)))
        print("Selecting 2 to {} prototypes.\n".format(len(labels) - 1))

        proto_dict = {}

        for k in range(2, len(labels)):
            # run prototypeSelection function
            prototypes = prototype_selection_destructive_maxdist(dm, k)
            proto_dict[k] = [labels[int(x)] for x in prototypes]

        with open(output[0], 'w') as outfile:
            dump(proto_dict, outfile, default_flow_style=False)

        with open(log[0], "w") as logfile:
            logfile.write("Running the run_prototypeSelection.py script.\n"
                          "The imported distance matrix has {0} elements.\n"
                          "Selecting 2 to "
                          "{1} prototypes.\n".format(df.shape[1],
                                                     len(labels) - 1))
