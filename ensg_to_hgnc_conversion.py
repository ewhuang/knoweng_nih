### Author: Edward Huang

from collections import OrderedDict
import numpy as np

### This script converts the new auc file for drug response in patients to the
### old format that the Mayo data used. Also converts the gene expression table
### to the old format. Lastly, converts the LINCS level 4 data to the old
### format.

if __name__ == '__main__':
    # Read in the cell line translation information.
    f = open('./data/cell_line.txt', 'r')
    # Keys are ccl_id's, values are the ccl_names.
    ccl_id_to_name = {}
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        # Skip a line if it doesn't have a translation.
        if len(line) == 1:
            continue
        ccl_name, ccl_id = line
        ccl_id_to_name[ccl_id] = ccl_name
    f.close()

    # Read in the drug translation information.
    f = open('./data/drug.txt', 'r')
    # Keys are master_cpd_id, values are broad_cpd_id.
    master_cpd_to_broad = {}
    for i, line in enumerate(f):
        if i == 0:
            continue
        master_cpd_id, broad_cpd_id = line.split()
        master_cpd_to_broad[master_cpd_id] = broad_cpd_id
    f.close()

    # Read in the raw drug response data, and construct the corresponding
    # dictionary.
    ccl_lst = []
    drug_resp_dct = OrderedDict({})
    f = open('./data/auc.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        experiment_id, auc, master_cpd_id, ccl_id = line.split()
        # Translate the drug id's and the ccl id's.
        drug_id = master_cpd_to_broad[master_cpd_id]
        # Insert values into drug response dictionary.
        if drug_id not in drug_resp_dct:
            drug_resp_dct[drug_id] = OrderedDict({})

        # Skip cancer cell lines not in our dictionary.
        if ccl_id not in ccl_id_to_name:
            continue
        ccl_name = ccl_id_to_name[ccl_id]
        if ccl_name not in ccl_lst:
            ccl_lst += [ccl_name]
        auc = float(auc)
        if ccl_name in drug_resp_dct[drug_id]:
            drug_resp_dct[drug_id][ccl_name] += [auc]
        else:
            drug_resp_dct[drug_id][ccl_name] = [auc]
    f.close()

    # Write out to our old format for drug response.
    out = open('./data/auc_hgnc.tsv', 'w')
    out.write('exposure\t' + '\t'.join(ccl_lst) + '\n')
    for drug in drug_resp_dct:
        out.write(drug)
        for ccl in ccl_lst:
            if ccl not in drug_resp_dct[drug]:
                out.write('\tNA')
            else:
                out.write('\t%f' % np.mean(drug_resp_dct[drug][ccl]))
        out.write('\n')
    out.close()

    f = open('./data/gene_expression.txt', 'r')
    out = open('./data/gene_expression_hgnc.tsv', 'w')
    # Gene expression cancer cell lines.
    for i, line in enumerate(f):
        if i == 0:
            out.write('gid\t')
            ge_ccl = line.split()[2:]
            # The set of cell lines in gene expression but not in drug response.
            ge_not_dr = set(ge_ccl).difference(ccl_lst)
            # Find indices of cell lines that should not be written out.
            ge_not_dr_indices = [i for i, e in enumerate(ge_ccl) if e not in ge_not_dr]
            ge_ccl = [ge_ccl[i] for i in ge_not_dr_indices]

            ge_indices = [ge_ccl.index(cl) for cl in ccl_lst]
            ge_ccl = [ge_ccl[i] for i in ge_indices]
            out.write('\t'.join(ge_ccl) + '\n')
            continue
        line = line.split()[1:]
        gene_id, cell_lines = line[0], line[1:]
        try:
            float(gene_id)
        except ValueError:
            # Remove cell lines that don't appear in the drug response.
            cell_lines = [cell_lines[i] for i in ge_indices]

            out.write(gene_id + '\t' + '\t'.join(cell_lines) + '\n')
    f.close()
    out.close()