from flask import Flask ,request, render_template
from decimal import Decimal
from timeit import default_timer as timer

needlemanWunsch = Flask(__name__)

@needlemanWunsch.route('/', methods=['GET','POST'])
def homepage():
        start = timer()
        # Open the file for reading 
        with open('seq1.txt','r') as file1 : lines1 = file1.readlines()
        with open('seq2.txt','r') as file2 : lines2 = file2.readlines()

        n = lines1[0]
        m = lines2[0]
        seq1 = lines1[1]
        seq2 = lines2[1]
        file1.close()
        file2.close()


        if request.method == "POST":
            match = request.form['match']
            misMatch = request.form['mis_match']
            gapVal = request.form['gapVal']
            aligned_seq1, aligned_seq2,scoreMatrix,seq1,seq2,score,match_score,mismatch_score,answerList = needleman_wunsch(seq1, seq2,int(n),int(m),int(match),int(misMatch),int(gapVal))
        else:
            aligned_seq1, aligned_seq2,scoreMatrix,seq1,seq2,score,match_score,mismatch_score,answerList = needleman_wunsch(seq1, seq2,int(n),int(m))
        arr = [aligned_seq1,aligned_seq2]
        end = timer()
        return render_template("home.html",mainarr = arr,score_matrix = scoreMatrix,seq1 = seq1, seq2 = seq2,
        score = score,match_score = match_score,mismatch_score = mismatch_score,
        tim = f"{((end - start) * 1000):.4f}", 
        answerList = answerList
        )

def needleman_wunsch(seq1, seq2,n,m, match_score=1, mismatch_score=-1, gap_score=1):
    # Initialize the score matrix
    score_matrix = [[0 for j in range(m+1)] for i in range(n+1)]

    # Initialize the first row and column of the score matrix
    for i in range(n+1):
        if i != 0:
            score_matrix[i][0] = score_matrix[i-1][0] + gap_score
        else:
            score_matrix[i][0] = 0
    for j in range(m+1):
        if j != 0:
            score_matrix[0][j] = score_matrix[0][j-1] + gap_score
        else:
            score_matrix[0][j] = 0

    # Fill in the rest of the score matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            state1 = score_matrix[i-1][j] + gap_score
            state2 = score_matrix[i][j-1] + gap_score
            score_matrix[i][j] = max(match, state1, state2)
    # Traceback and build the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    j = m
    i = n
    score = 0
    answersList = [[1,2]]
    while i > 0 or j > 0:
        score_diag = score_matrix[i-1][j-1] if i > 0 and j > 0 else Decimal('-inf')
        score_up = score_matrix[i-1][j] if i > 0 else Decimal('-inf')
        indexUp = [i-1,j]
        score_left = score_matrix[i][j-1] if j > 0 else Decimal('-inf')
        indexLeft = [i,j-1]
        maxVal = max(score_diag,score_left,score_up)
        if maxVal == score_diag or (seq1[i-1] == seq2[j-1] and len(seq1) == len(seq2)):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            answersList.append([i,j])
            i -= 1
            j -= 1
        elif (score_left == score_up and maxVal != score_diag):
            # print(seq1[i - 1])
            left = score_left + (match_score if seq1[indexLeft[0]-1] == seq2[indexLeft[1]-1] else mismatch_score)
            up = score_up + (match_score if seq1[indexUp[0]-1] == seq2[indexUp[1]-1] else mismatch_score)
            # print(left)
            if (left > up):
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                answersList.append([i,j])
                j -= 1
            else:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                answersList.append([i,j])
                i -= 1
        elif maxVal == score_left :
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            answersList.append([i,j])
            j -= 1
        elif maxVal == score_up:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            answersList.append([i,j])
            i -= 1
    
    # Calculate the score
    for i in range(0,len(aligned_seq1),1):
        if aligned_seq1[i] == aligned_seq2[i]:
            score  = score + match_score
        elif aligned_seq1[i] == '-' or  aligned_seq2[i] == '-':
            score  = score + gap_score
        else:
            score = score + mismatch_score

    return (aligned_seq1, aligned_seq2,score_matrix,seq1,seq2,score,match_score,mismatch_score,answersList)




if __name__ == "__main__":
    
    needlemanWunsch.run(debug=True,port=9000)