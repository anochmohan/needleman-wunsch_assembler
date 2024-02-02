#!/usr/bin/env python3

import numpy as np
import sys
def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    '''Create Needleman-Wunsch algorithm
        seq_a = first seqeuence
        seq_b = second sequence
        match = score given for matches
        mismatch = score given for mismatches
        gap = score given for mismatches'''
    
    ### CONSTRUCT MATRIX USING seq_a and seq_b
    matrix = np.zeros((len(seq_a)+1, len(seq_b)+1))
    #print(matrix)

    # Fill top row and left column as gaps
    for i in range(len(seq_a)+1):
        matrix[i][0] = i*gap
    for j in range(len(seq_b)+1):
        matrix[0][j] = j*gap

    #print(matrix)

    # ABOVE WE CREATED A  MATRIX[i,j], where the first row and first column are cumulative gaps (0,-1,-2,...)
    # THISS MEANS THAT THE ACTUAL VALUES CORRESPPONDING TO THE NUCLEOTIDES ARE ALWAYS GOING TO BE OFF BY -1 

    ### FILL THE MATRIX
    # NOW WE WILL START FILLING THE MATRIX
    # THERE ARE 3 PARTS TO FILLING THE MATRIX:
        ### 1) IF A CELL IS (i-1, j-1) - this means the row and column had a nucleotide. (match bonus / mismatch penalty)
        ### 2) IF A CELL IS (i-1, j)   - this means the row matches, but column has insertion. (gap penalty)
        ### 3) IF A CELL IS (i, j-1)   - this means the column matches, but row has deletion. (gap penalty)

    # Matrix to check score
    matrix_with_score = np.zeros((len(seq_a), len(seq_b)))

    # Add scores based on match/mismatch to matrix_with_score
    for i in range(len(seq_a)):
        for j in range(len(seq_b)):
            if seq_a[i] == seq_b[j]:
                matrix_with_score[i][j] = match
            else:
                matrix_with_score[i][j] = mismatch            
    #print(matrix_with_score)

    # FILL THE MATRIX WITH SCORES
    for i in range(1, len(seq_a)+1):
        for j in range(1, len(seq_b)+1):
            nucleotide_matched = matrix[i-1][j-1]+matrix_with_score[i-1][j-1] # 1) insert match/mismatch score
            deletion = matrix[i-1][j] + gap # 2) insert gap score for insersions
            insertion = matrix[i][j-1] + gap # 3) insert gap score for deletions
            matrix[i][j] = max(nucleotide_matched, deletion, insertion)
    #print(matrix)

    # Trace the best pathh back 
    def traceback(matrices):    
        new_seq_a = []
        new_seq_b = []
        num_rows = len(seq_a)
        num_cols = len(seq_b)

        current_row, current_col = num_rows, num_cols

        while current_row > 0 and current_col > 0:
            current_score = matrix[current_row][current_col]
            up_score = matrix[current_row-1][current_col]
            #left_score = matrix[current_row][current_col-1]
            diagonal_score = matrix[current_row-1][current_col-1]
            if diagonal_score == current_score - match and seq_a[current_row-1] == seq_b[current_col-1]:
                new_seq_a.append(seq_a[current_row-1])
                new_seq_b.append(seq_b[current_col-1])
                current_row -= 1
                current_col -= 1
            elif diagonal_score == current_score - mismatch and seq_a[current_row-1] != seq_b[current_col-1]:
                new_seq_a.append(seq_a[current_row-1])
                new_seq_b.append(seq_b[current_col-1])
                current_row -= 1
                current_col -= 1
            elif up_score == current_score - gap:
                new_seq_b.append('-')
                new_seq_a.append(seq_a[current_row-1])
                current_row -= 1
            else:
                new_seq_a.append('-')
                new_seq_b.append(seq_b[current_col-1])
                current_col -= 1
            score = matrix[-1][-1]
        # Reverse the lists to get the correct order
        new_seq_a.reverse()
        new_seq_b.reverse()
        new_seq_a = "".join(new_seq_a)
        new_seq_b = "".join(new_seq_b)

        return new_seq_a, new_seq_b, int(score)
    
    final_seq_a, final_seq_b, final_score = traceback(matrix)

    return  ((final_seq_a, final_seq_b), final_score)