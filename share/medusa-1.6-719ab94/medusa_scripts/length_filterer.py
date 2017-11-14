from sys import argv,exit

def filtering(inp,out="filtered.fna"):
	from Bio import SeqIO
	handle = open(inp, "rU")
	goods=[r for r in SeqIO.parse(handle,'fasta') if len(r) >= 1000]
	sequences = goods
	output_handle = open(out, "w")
	SeqIO.write(sequences, output_handle, "fasta")
	output_handle.close()

if __name__ == '__main__':
	if len(argv) == 2: inp=argv[1]
	elif len(argv) == 3:
		inp,out=argv[1:]
	else:
		print 'wrong usage'
		exit()
	try:filtering(inp,out)
	except:filtering(inp)

