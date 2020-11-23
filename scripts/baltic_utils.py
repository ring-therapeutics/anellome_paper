import importlib.machinery
import importlib.util
loader = importlib.machinery.SourceFileLoader('baltic','/Users/evogytis/Documents/baltic/baltic/baltic.py')
spec = importlib.util.spec_from_loader(loader.name, loader)
bt = importlib.util.module_from_spec(spec)
loader.exec_module(bt)

from Bio import SeqIO,Seq ## Biopython
import os,glob,re

def reconstruct_ancestral_states(cfml_path,aln_path,tree_path,kappa,out_stem,include_seq=True,traits=None):
    """
    Given the path to ClonalFrameML, a fasta alignment, a newick tree, a kappa and an output stem runs ClonalFrameML.
    Produces intermediate renamed files and ClonalFrameML files and returns a baltic tree with reconstructed sequences on branches.
    """
    nonnucleotide=re.compile('[^ACTG\-Nactgn]') ## recognises nucleotides that aren't A,C,T,G,N or - (ambiguous codes)

    rename=[]
    newNames={}
    oldNames={}
    for i,seq in enumerate(SeqIO.parse(aln_path,format='fasta')): ## iterate over sequences in fasta
        newNames[seq.id]=str(i) ## new name is just index
        oldNames[str(i)]=seq.id ## remember old name
        seq.id=str(i) ## assign index to sequence instead of name (ClonalFrameML doesn't like symbols in sequence names)
        seq.description=''

        S=str(seq.seq) ## turn sequence to string

        while nonnucleotide.search(S): ## search and replace ambiguous nucleotide codes
            S=nonnucleotide.sub('N',S) ## replace all ambiguous sites with Ns

        seq.seq=Seq.Seq(S) ## assign new sequence, free of ambiguous nucleotides
        rename.append(seq)

    fasta=aln_path.replace('.fasta','_renamed.fasta') ## new fasta name
    out=open(fasta,'w') ## start empty file
    for seq in rename: ## iterate over sequences loaded previously
        out.write('>%s\n%s\n'%(seq.id,seq.seq)) ## write out a renamed fasta entry into new file
    out.close() ## close file

    ll=bt.loadNewick(tree_path,absoluteTime=False) ## load tree file
    ll.treeStats() ## report some stats

    for k in ll.getExternal(): ## iterate over tips in the tree
        k.name=newNames[k.name] ## assign new names to tips

    newick=tree_path.replace('.newick','_renamed.newick') ## create new file name
    out=open(newick,'w') ## open empty tree file
    out.write('%s\n'%(ll.toString().replace("'",""))) ## write out the newick string to file
    out.close() ## close tree file

    base_path=os.path.dirname(newick)

    print('\n%s %s %s %s -imputation_only true -kappa %s\n'%(cfml_path,newick,fasta,out_stem,kappa))
    os.system('cd %s; %s %s %s %s -imputation_only true -kappa %s'%(base_path,cfml_path,newick,fasta,out_stem,kappa)) ## call ClonalFrameML on fasta and newick files with renamed tips

    cross_ref=open(os.path.join(base_path,out_stem+'.position_cross_reference.txt'),'r').readline() ## load position cross reference file from reconstruction
    cross_ref={i:x for i,x in enumerate(map(int,cross_ref.split(',')))} ## index position references

    reconstructed_SNPs=open(os.path.join(base_path,out_stem+'.ML_sequence.fasta'),'r') ## load reconstructed SNPs
    original_name_aln=open(os.path.join(base_path,out_stem+'.reconstructed.fasta'),'w')

    SNPs={} ## will contain SNPs for each tip

    for line in reconstructed_SNPs: ## iterate over file
        l=line.strip('\n')
        if '>' in l: ## fasta header
            name=l[1:] ## name is on this line
        else:
            SNPs[name]=l ## remember sequence

    base_seq=[] ## contains invariant sites
    for line in open(fasta,'r'): ## iterate over original sequences
        l=line.strip('\n')
        if '>' in l: ## fasta header
            name=l[1:] ## name is on this line
        else: ## not fasta header, therefore sequence
            seq=list(l) ## turn sequence into a list of nucleotides
            for i in cross_ref: ## iterate over cross-reference
                if cross_ref[i]>0: ## variable site
                    seq[i]=SNPs[name][cross_ref[i]-1] ## fetch appropriate site from SNPs
                    base_seq.append('N') ## invariant sequence doesn't have anything here
                else: ## invariant site
                    base_seq.append(seq[i]) ## invariant sequence has the same nucleotide (non-variable site)
            seq=''.join(seq) ## concatenate sequence

            if 'NODE_' not in name: name=oldNames[name] ## tip, rename
            original_name_aln.write('>%s\n%s\n'%(name,seq))

    ll=bt.loadNewick(os.path.join(base_path,out_stem+'.labelled_tree.newick'),absoluteTime=False) ## load labelled tree
    for k in ll.getInternal(): ## iterate over nodes
        original_name_aln.write('>%s\n%s\n'%(k.traits['label'],''.join([base_seq[i] if cross_ref[i]==0 else SNPs[k.traits['label']][cross_ref[i]-1] for i in cross_ref]))) ## attach sequence to branch
    for k in ll.getExternal():
        k.name=oldNames[k.name]

    backname_tree=open(os.path.join(base_path,out_stem+'.labelled_tree.newick'),'w')
    backname_tree.write(ll.toString())
    backname_tree.close()



    original_name_aln.close() ## close alignment
    # original_name_aln=open(os.path.join(base_path,out_stem+'.ML_sequence.original_name.fasta'),'w')


    # base_seq=[] ## contains invariant sites
    # for line in open(fasta,'r'): ## iterate over original sequences
    #     l=line.strip('\n')
    #     if '>' in l: ## fasta header
    #         name=l[1:] ## get name
    #         original_name_aln.write('>%s\n'%(oldNames[int(l[1:])])) ## rename sequence in new file
    #     else: ## sequence
    #         seq=list(line) ## turn sequence into a list of nucleotides
    #         for i in cross_ref: ## iterate over cross-reference
    #             if cross_ref[i]>0: ## variable site
    #                 seq[i]=SNPs[name][cross_ref[i]-1] ## fetch appropriate site from SNPs
    #                 base_seq.append('N') ## invariant sequence doesn't have anything here
    #             else:
    #                 base_seq.append(seq[i]) ## invariant sequence has the same nucleotide (non-variable site)
    #         seq=''.join(seq) ## concatenate sequence
    #         original_name_aln.write(seq) ## write alignment to file
    #         reconstructed_seqs[name]=seq ## remember reconstructed sequence

    # original_name_aln.close() ## close alignment


    # ll=bt.loadNewick(os.path.join(base_path,out_stem+'.labelled_tree.newick'),absoluteTime=False) ## load labelled tree
    # original_name_tree=open(os.path.join(base_path,out_stem+'.annotated.tree'),'w') ## open new file
    #
    # for k in ll.getExternal(): ## iterate over tips
    #     k.traits['sequence']=reconstructed_seqs[k.name].strip('\n') ## add sequence
    #     k.name=oldNames[int(k.name)] ## rename
    # for k in ll.getInternal(): ## iterate over nodes
    #     k.traits['sequence']=''.join([base_seq[i] if cross_ref[i]==0 else SNPs[k.traits['label']][cross_ref[i]-1] for i in cross_ref]) ## attach sequence to branch
    #
    # include_traits=['label']
    # if traits!=None:
    #     include_traits+=traits
    # if include_seq==True: include_traits.append('sequence')
    #
    # original_name_tree.write('%s'%(ll.toString(traits=include_traits,nexus=True))) ## write tree out with node labels, but no sequences
    # original_name_tree.close() ## close tree
    #
    # mutation_file=open(os.path.join(base_path,out_stem+'.mutations.tsv'),'w') ## open new file
    # for k in ll.Objects: ## iterate over branches
    #     name=k.name if k.branchType=='leaf' else k.traits['label']
    #     line=[name]
    #     if 'sequence' in k.parent.traits: ## parent isn't root
    #         seq=k.traits['sequence'] ## get current sequence
    #         par_seq=k.parent.traits['sequence'] ## get parent sequence
    #         par_name=k.parent.traits['label']
    #         mutations=[]
    #         for site,pair in enumerate(zip(seq,par_seq)): ## iterate over every pair of sites
    #             A,B=pair ## get child and parent states
    #             if A!=B: ## they don't match - there's been a mutation
    #                 if 'mutations' not in k.traits: k.traits['mutations']=[] ## start mutation list
    #                 k.traits['mutations'].append('%s%d%s'%(A,site+1,B)) ## assign mutation (site indexing is 1-based, not 0-based like python)
    #                 mutations.append('%s%d%s'%(A,site+1,B))
    #
    #         line.append(','.join(mutations))
    #         line.append(par_name)
    #     mutation_file.write('%s\n'%('\t'.join(line)))
    # mutation_file.close()
    # return ll

def annotate_tree(tree,tree_file,seq_file):
    seqs={}
    for seq in SeqIO.parse(seq_file,'fasta'):
        seqs[seq.id]=seq.seq

    ll=bt.loadNewick(os.path.join(tree_file),absoluteTime=False) ## load labelled tree
    relationships={}
    name=lambda k: k.name if k.branchType=='leaf' else k.traits['label']

    for k in ll.Objects:
        if 'label' in k.parent.traits: relationships[name(k)]=name(k.parent)

    for k in tree.Objects:
        if 'label' in k.traits: del k.traits['label'] ## reset labels in case the tree somehow had them before

    while any(['label' not in k.traits for k in tree.getInternal()]): ## annotate internal nodes
        for k in tree.Objects:
            if ((k.branchType=='node' and 'label' in k.traits) or k.branchType=='leaf') and name(k) in relationships: ## (node has label or we have tip) and its relationship is known
                k.parent.traits['label']=relationships[name(k)]

    for k in tree.Objects:
        k.traits['seq']=seqs[name(k)] ## annotate sequence on node

    for k in tree.Objects: ## iterate over branches
        if 'seq' in k.parent.traits: ## parent isn't root
            seq=k.traits['seq'] ## get current sequence
            par_seq=k.parent.traits['seq'] ## get parent sequence
            for site,pair in enumerate(zip(seq,par_seq)): ## iterate over every pair of sites
                A,B=pair ## get child and parent states
                if A!=B: ## they don't match - there's been a mutation
                    if 'muts' not in k.traits: k.traits['muts']=[] ## start mutation list
                    k.traits['muts'].append('%s%d%s'%(A,site+1,B)) ## assign mutation (site indexing is 1-based, not 0-based like python)

    # return tree
