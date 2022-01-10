from Bio import SeqIO
def getAverage_getConsenesus():
    Lines=[]#initial array for storing sequences
    with open("averages & consenesus.txt","w") as f:#make a new txt file and open it
        for seq_record in SeqIO.parse("gisaid_hcov-19_2021_12_31_12.fas", "fasta"):#get the data from fasta file
                #initial variables------------------------------
                vertical_Seq=[]
                vertical_sequences_list=[]
                global new_Sequence
                new_Sequence=[]
                element_counts=[]
                #--------------------------------------------------
                Lines.append(str(seq_record.seq))
                #debuging part-------------------------
                mysequence=str(seq_record.seq)
                #-------------------------------------
                #calculating the averages of delta sequences----------------------------
                avg_a=100*(str(seq_record.seq).count('A')/len(str(seq_record.seq)))
                avg_g=100*(str(seq_record.seq).count('G')/len(str(seq_record.seq)))
                avg_c=100*(str(seq_record.seq).count('C')/len(str(seq_record.seq)))
                avg_t=100*(str(seq_record.seq).count('T')/len(str(seq_record.seq)))
                avg_cg=avg_c+avg_g
                #printing them to the txt file---------------------
                f.write("Sequence averages values: " + "\n")
                f.write("average of (A)= "+str(avg_a)+"%" + "\n")
                f.write("average of (G)= "+str(avg_g )+"%"+ "\n")
                f.write("average of (C)= "+str(avg_c)+"%" + "\n")
                f.write("average of (T)= "+str(avg_t )+"%"+ "\n")
                f.write("average of (CG)= "+str(avg_cg)+"%" + "\n")
                f.write("---------------------------------------------------------------------------------"+"\n")

        for i in range(len(Lines[0])-1):#storing vertical colums
            for line in Lines:
                vertical_Seq.append(line[i])
            vertical_sequences_list.append(vertical_Seq)
            vertical_Seq=[]
        
        for i in range(len(Lines[0])-1):#storing the count of every element of colum and store the maximum one to consensus seq
            for j in range (10):
                elementCount=vertical_sequences_list[i].count(vertical_sequences_list[i][j])
                element_counts.append(elementCount)
            index=element_counts.index(max(element_counts))
            new_Sequence.append(vertical_sequences_list[i][index])
            element_counts=[]
        #calculating the averages of consensus sequences----------
        new_Sequence=''.join(new_Sequence)
        cons_avg_a=100*(str(new_Sequence).count('A')/len(new_Sequence))
        cons_avg_g=100*(str(new_Sequence).count('G')/len(new_Sequence))
        cons_avg_c=100*(str(new_Sequence).count('C')/len(new_Sequence))
        cons_avg_t=100*(str(new_Sequence).count('T')/len(new_Sequence))
        cons_avg_cg=cons_avg_c+cons_avg_g
        #printing the consensus sequence and its avg values to the file 
        f.write("consensus Sequence: " + "\n")
        f.write(str(new_Sequence) + "\n")
        f.write("consensus Sequence averages values: "  +"\n")
        f.write("average of (A)= "+str(cons_avg_a)+"%" + "\n")
        f.write("average of (G)= "+str(cons_avg_g)+"%" + "\n")
        f.write("average of (C)= "+str(cons_avg_c)+"%"+ "\n")
        f.write("average of (T)= "+str(cons_avg_t )+"%"+ "\n")
        f.write("average of (CG)= "+str(cons_avg_cg)+"%" + "\n")
        
                
def comparing_sequences():
    Lines=[]
    print(len(new_Sequence))
    with open("comparing results.txt","w") as f:#make a new txt file and open it
        f.write("the indices of different bases of omecron are: "+ "\n")
        for seq_record in SeqIO.parse("gisaid_hcov_omicron-19_2021_12_31_12.fas", "fasta"):#get the data 
            Lines.append(str(seq_record.seq))
            mysequence=str(seq_record.seq)
        vertical_Seq=[]
        vertical_sequences_list=[]
        element_counts=[]
        for i in range(len(Lines[0])-1):#storing the colums again but now of omecron
            for line in Lines:
                vertical_Seq.append(line[i])
            vertical_sequences_list.append(vertical_Seq)
            vertical_Seq=[]
        indices=[]
        for i in range(len(Lines[0])-1): #storing the indices of elements that repeated at least 7 times
            for j in range (10):
                elementCount=vertical_sequences_list[i].count(vertical_sequences_list[i][j])
                if elementCount>6:
                    index=i
                    indices.append(index)
        indices = list( dict.fromkeys(indices) )      
        for index in indices:#compaing the elements of these indices----prints when there is disimilar
            if Lines[0][index]!=new_Sequence[index]:
                f.write("index ("+str(index)+")"+" was "+new_Sequence[index]+" became "+mysequence[index]+" - ")
                f.write("\n")
                f.write("*************************************************************************"+"\n")
            
            
getAverage_getConsenesus()
comparing_sequences()
