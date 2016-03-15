

import math
import operator

# def min_p_exp(p_val_lst):
#     return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))


# drug_list = []
# f = open('./results/best_drugs_ppi_lincs.txt', 'r')
# for line in f:
#     line = line.strip()
#     drug_list += [line]
# f.close()


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

# # print 'Extracting level 4 LINCS top pathways...'
# f = open('./results/top_pathways_lincs_Aft_%s_hgnc.txt' % 3, 'r')
# lincs_dct = {}
# for i, line in enumerate(f):
#     # Skip first two lines, they're just results summaries.
#     if i < 2:
#         continue

#     # Disregard x, y, z.
#     drug, cell_line, path, score, x, y, z = line.strip().split('\t')
#     if drug not in drug_list:
#         continue
#     score = float(score)

#     if (drug, path) not in lincs_dct:
#         lincs_dct[(drug, path)] = [score]
#     else:
#         lincs_dct[(drug, path)] += [score]
# f.close()

# # Aggregate the p-values across cell lines for each drug-pathway pair.          
# for (drug, path) in lincs_dct:
#     assert drug in drug_list
#     p_val_lst = lincs_dct[(drug, path)]
#     # Change the sum_log_p_values() function for other aggregation
#     # functions.
#     drug_path_p_val = min_p_exp(p_val_lst)
#     lincs_dct[(drug, path)] = drug_path_p_val

ppi_dct = {}
f = open('./results/embedding/genetic_top_pathways_100_0.8.U_top_250.txt', 'r')
for i, line in enumerate(f):
    if i < 2:
        continue
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