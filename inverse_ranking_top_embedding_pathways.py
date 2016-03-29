
import file_operations
import math
import operator
import sys

def get_inverse_rankings(subfolder, filename):   
    brd_drug_to_name_dct = file_operations.get_brd_drug_to_name_dct()

    print 'Taking out top 20 global pathways'
    top_global_pathways = file_operations.get_top_global_pathways()[:20]

    drug_path_dct = {}
    f = open(subfolder + filename, 'r')
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug, path, score = line
        if path in top_global_pathways:
            continue
        drug_path_dct[(drug, path)] = float(score)
    f.close()

    # Sort ppi dictionary by value.
    ranked_drug_path_dct = {}
    drug_path_dct = sorted(drug_path_dct.items(), key=operator.itemgetter(1),
        reverse=True)
    for i, ((drug, path), score) in enumerate(drug_path_dct):
        inverse_ranking = (i + 1) / float(len(drug_path_dct))
        ranked_drug_path_dct[(drug, path)] = inverse_ranking

    ranked_drug_path_dct = sorted(ranked_drug_path_dct.items(),
        key=operator.itemgetter(1), reverse=True)

    out = open(subfolder + 'inverse_' + filename, 'w')
    out.write('drug\tpath\tinverse_rank\n')
    for (drug, path), score in ranked_drug_path_dct:
        out.write('%s\t%s\t%g\n' % (brd_drug_to_name_dct[drug], path, score))
    out.close()

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s method dimension U_type' % sys.argv[0]
        exit()
    method = sys.argv[1]
    dimension = sys.argv[2]
    U_type = sys.argv[3]
    assert method in ['ppi', 'genetic', 'literome', 'sequence']
    assert dimension in ['50', '100', '500', '1000', '1500', '2000']
    assert U_type in ['U', 'US']
    filename = '%s_top_pathways_%s_0.8.%s_top_250.txt' % (method, dimension,
        U_type)
    subfolder = './results/embedding/'
    get_inverse_rankings(subfolder, filename)

if __name__ == '__main__':
    main()