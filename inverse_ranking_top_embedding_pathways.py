

import math
import operator

drug_to_name_dct = {}
f = open('./data/brd_to_name.txt', 'r')
for i, line in enumerate(f):
    # Skip header line.
    if i == 0 or line == '\n':
        continue
    line = line.strip().split('\t')
    assert len(line) == 3
    drug_name = line[0]
    drug_id = line[2]
    assert drug_id not in drug_to_name_dct
    drug_to_name_dct[drug_id] = drug_name
f.close()

ppi_dct = {}
f = open('./results/top_pathways_genetic_subtract_superdrug.txt', 'r')
for i, line in enumerate(f):
    # if i < 2:
    #     continue
    line = line.strip().split('\t')
    assert len(line) == 3
    drug, path, score = line
    # if drug not in drug_list:
    #     continue
    score = float(score)
    ppi_dct[(drug, path)] = score
f.close()

# Sort ppi dictionary by value.
ranked_ppi_dct = {}
ppi_dct = sorted(ppi_dct.items(), key=operator.itemgetter(1), reverse=True)
num_pairs = float(len(ppi_dct))
for i, ((drug, path), score) in enumerate(ppi_dct):
    new_score = (i + 1) / num_pairs
    ranked_ppi_dct[(drug, path)] = new_score

out = open('./results/genetic_inverse_rank.txt', 'w')
# out.write('drug\tpath\tlincs_p_value\tppi_inverse_rank\n')
out.write('drug\tpath\tgenetic_inverse_rank\n')
# lincs_dct = sorted(lincs_dct.items(), key=operator.itemgetter(1))
ranked_ppi_dct = sorted(ranked_ppi_dct.items(), key=operator.itemgetter(1),
    reverse=True)
# for (drug, path), score in lincs_dct:
    # assert (drug, path) in ranked_ppi_dct
for (drug, path), score in ranked_ppi_dct:
    # out.write('%s\t%s\t%g\t%g\n' % (drug_to_name_dct[drug], path, score, 
    #     ranked_ppi_dct[(drug, path)]))
    out.write('%s\t%s\t%g\n' % (drug_to_name_dct[drug], path, score))


out.close()