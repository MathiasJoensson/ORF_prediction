#create GUI window for DNA Sequence Analyzer

from cgitb import text
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from DNAToolkit import *

#Create window
window = Tk()
window.title("DNA Sequence Analyzer")
window.geometry("500x500")
window.configure(bg="white")

#Create label
label = Label(window, text="DNA Sequence Analyzer", font=("Arial Bold", 20), bg="white")
label.pack(pady=10)

#Create a textbox
textbox = Text(window, height=10, width=40, font = ('Courier', 12))
textbox.pack(pady=10)

#Create a label for the result
result_label = Label(window, text="", font=("Arial", 12), bg="white")
result_label.pack(pady=10)

#Create a button to browse files
def browseFiles():
    filename = filedialog.askopenfilename(initialdir = "/", title = "Select a File", filetypes = (("Text files", "*.txt*"), ("all files", "*.*")))
    if filename:
        with open(filename, "r") as f:
            text.delete(1.0, END)
            text.insert(1.0, f.read())

button = Button(window, text="Browse Files", command=browseFiles)
button.pack(pady=10)

#Create a function to display the results in a scrollable container
def displayResults(results):
    #Create a new window to display the results
    newWindow = Toplevel(window)
    newWindow.title("Results")
    newWindow.geometry("500x500")
    newWindow.configure(bg="white")

    #Create a label for the result
    result_label = Label(newWindow, text="Results", font=("Arial", 12), bg="white")
    result_label.pack(pady=10)

    #Create a scrollable container
    scrollbar = Scrollbar(newWindow)
    scrollbar.pack(side=RIGHT, fill=Y)

    #Create a container to display the results
    result_box = Text(newWindow, height=10, width=40, font = ('Courier', 12), yscrollcommand=scrollbar.set)
    result_box.pack(pady=10)

    #Insert the results into the container
    for result in results:
        result_box.insert(END, result)
        result_box.insert(END, "\n")
    
    #Attach the scrollbar to the container
    scrollbar.config(command=result_box.yview)




#Create labels for the results
freq_label = Label(window, text="", font=("Arial", 12), bg="white")
gc_label = Label(window, text="", font=("Arial", 12), bg="white")
rna_label = Label(window, text="", font=("Arial", 12), bg="white")
aa_label = Label(window, text="", font=("Arial", 12), bg="white")
frame_label = Label(window, text="", font=("Arial", 12), bg="white")
proteins_label = Label(window, text="", font=("Arial", 12), bg="white")
longest_label = Label(window, text="", font=("Arial", 12), bg="white")

#Create a button to submit sequence
def submit():
    seq = textbox.get(1.0, END).strip()
    print(seq)
    if seq:
        #Call function to validate sequence
        valid_seq = validateSeq(seq)
        if valid_seq:
            result_label.config(text="The sequence is valid")
            #Call other functions to analyze the sequence
            #Count nucleotide frequency
            freq = countNucFrequency(valid_seq)
            freq_label.config(text=f'Nucleotide frequency: {freq}')

            #TRANSCRIPTION
            #DNA -> RNA Transcription
            rna = transcription(valid_seq)
            rna_label.config(text=f'DNA -> RNA transcription: {rna}')

            #DNA Reverse compliment
            #rev = reverse_compliment(valid_seq)

            #GC content calculation
            gc = gc_content(valid_seq)
            gc_label.config(text=f'GC content: {gc}%')

            #TRANSLATION
            #translate DNA sequence into AA sequence, init_pos indicates the reading frame
            aa = translate_seq(valid_seq, 0)
            aa_label.config(text=f'Amino acid seq from DNA: {aa}')

            #Generate 6 reading frames
            frame = gen_reading_frames(valid_seq)
            frame_label.config(text=f'Reading frames: {frame}')

            #protein search within the reading frames
            proteins = all_proteins_from_orfs(valid_seq, 0, 0, True)
            proteins_label.config(text=f'Proteins: {proteins}')

            #Longest protein sequence
            longest = longest_protein_from_orfs(valid_seq, 0, 0)
            longest_label.config(text=f'Longest protein: {longest}')
            
            #Display results in a new window
            results = [
                f'\nDNA sequence: {valid_seq}\n',
                f'Sequence Length: {len(valid_seq)}\n',
                f'Nucleotide frequency: {freq}\n',
                f'DNA -> RNA transcription: {rna}\n',
                f'GC content: {gc}%\n',
                f'Amino acid seq from DNA: {aa}\n',
                f'Reading frames: {frame}\n',
                f"Proteins: {proteins}\n",
                f'Longest protein: {longest}\n'
            ] 
            displayResults(results)
            print(f'{proteins}')

        else:
            result_label.config(text="The sequence is not valid")
    else:
        #Show error message here
        result_label.config(text="Please enter or upload a sequence")

button = Button(window, text="Submit", command=submit, font = ('Arial', 12), bg="#d3d3d3")
button.pack(pady=10)


#Run the main loop
window.mainloop()