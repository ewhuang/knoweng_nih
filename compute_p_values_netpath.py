### Author: Edward Huang

import operator
import time

### This script compares the 10k runs of random embedding to the proper run of
### embedding, considering cosine scores only. The number of cosine scores in
### the random run that are higher than the proper run for each pathway is the
### pseudo-p-value for that drug-pathway pair.

def read_random_top_pathways():
    '''
    Returns a dictionary.
    Key: str -> pathway name.
    Value: [floats] -> cosine scores from random embedding experiment..
    '''
    random_pathway_dct = {}
    subfolder = './results/random_embedding/'
    f = open('%srandom_ppi_top_pathways_50_0.8.U_top_250.txt' % subfolder, 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        line = line.strip().split('\t')
        pathway, cosine_scores = line[0], map(float, line[1:])

        assert pathway not in random_pathway_dct
        random_pathway_dct[pathway] = cosine_scores
    f.close()
    return random_pathway_dct

def get_embedding_path_dct():
    '''
    Given the correlation method path dictionary, find the top pathways for
    each drug as computed by embedding, where the number of top pathways for
    each drug is the same as the number as computed by correlation.
    '''
    embedding_path_dct = {}

    f = open('./results/random_embedding/ppi_top_pathways_50_0.8.U_top_250_just_cosine.txt', 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        drug, path, score = line.strip().split('\t')
        if score == '[]':
            continue

        assert (drug, path) not in embedding_path_dct
        embedding_path_dct[(drug, path)] = float(score)
    f.close()
    return embedding_path_dct

def get_p_value_dct(random_pathway_dct, embedding_path_dct):
    '''
    Returns a dictionary.
    Key: (str, str) -> (drug, pathway)
    Value: int -> number of random cosine scores out of 10k better than the 
            netpath score for the same key.
    '''
    netpath_p_value_dct = {}
    for (drug, path) in embedding_path_dct:
        random_cosine_scores = random_pathway_dct[path]
        netpath_score = embedding_path_dct[(drug, path)]
        num_better = sum(i > netpath_score for i in random_cosine_scores)
        
        assert (drug, path) not in netpath_p_value_dct
        netpath_p_value_dct[(drug, path)] = num_better / float(len(
            random_cosine_scores))
    return netpath_p_value_dct 

def write_out_netpath_p_values(netpath_p_value_dct):
    '''
    Takes the NetPath p-values and writes them out to file.
    '''
    netpath_p_value_dct = sorted(netpath_p_value_dct.items(),
        key=operator.itemgetter(1))
    out = open('./results/random_embedding/netpath_p_values.txt', 'w')
    out.write('drug\tpath\tpseudo_p_value\n')
    for (drug, path), pseudo_p_val in netpath_p_value_dct:
        out.write('%s\t%s\t%f\n' % (drug, path, pseudo_p_val))
    out.close()

def main():
    random_pathway_dct = read_random_top_pathways()
    embedding_path_dct = get_embedding_path_dct()
    netpath_p_value_dct = get_p_value_dct(random_pathway_dct,
        embedding_path_dct)
    write_out_netpath_p_values(netpath_p_value_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))