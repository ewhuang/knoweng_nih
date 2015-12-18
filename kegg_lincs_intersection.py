### Author: Edward Huang

### This file finds the number of unique genes in the NCI pathways and the
### number of unique genes in the KEGG pathways. It prints out these numbers
### and also the size of their intersection.

if __name__ == '__main__':
    print 'Extracting NCI pathways...'
    path_file = open('./data/nci_pathway.txt', 'r')
    nci_pathnames = []
    nci_path_dct = {}
    nci_path_genes = set([])
    for line in path_file:
        line = line.split('\t')
        path_name, path_gene = line[0], line[1][:-2]
        nci_path_genes.add(path_gene)
        if path_name not in nci_path_dct:
            nci_path_dct[path_name] = [path_gene]
        else:
            nci_path_dct[path_name] += [path_gene]
        if path_name not in nci_pathnames:
            nci_pathnames += [path_name]
    path_file.close()


    print 'Extracting NCI pathways...'
    path_file = open('./data/kegg_pathway.txt', 'r')
    kegg_pathnames = []
    kegg_path_dct = {}
    kegg_path_genes = set([])
    for line in path_file:
        line = line.split('\t')
        path_name, path_gene = line[0], line[1]
        kegg_path_genes.add(path_gene)
        if path_name not in kegg_path_dct:
            kegg_path_dct[path_name] = [path_gene]
        else:
            kegg_path_dct[path_name] += [path_gene]
        if path_name not in kegg_pathnames:
            kegg_pathnames += [path_name]
    path_file.close()

    print len(kegg_path_genes), len(nci_path_genes), len(nci_path_genes.intersection(kegg_path_genes))