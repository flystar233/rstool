import csv
import re
import click
from  collections import defaultdict
import matplotlib.pyplot as plt

pattern = {'GT':'K','AC':'M','AG':'R','CG':'S','CT':'Y','AT':'W','TG':'K','CA':'M','GA':'R','GC':'S','TC':'Y','TA':'W'}
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('##')[0].strip()
        if raw: yield raw

@click.command("vcf_stat", help="stat vcf file")
@click.option(
    "-f", "--vcf", type=click.Path(exists=True), help="vcf file"
    )
@click.option(
    "-s", "--select", default="a", type=str, help="The number of alleles is one(a) or two(other)", show_default=True)
@click.option(
    "-v", "--outfile1", default="out.vcf",type=str, help="streamlined vcf file",show_default=True)
@click.option(
    "-g", "--outfile2", default="out.geno", type=str, help="genotype file ", show_default=True)
def vcf_stat(vcf,outfile1,outfile2,select):
	try:
		with open(vcf) as IN,open(outfile1,'w') as OUT1,open(outfile2,'w') as OUT2:
			f_csv = csv.reader(decomment(IN),delimiter='\t') # skip annoation(##)
			header = next(f_csv) # skip annoation(#)
			CHR = header[0]
			POS = header[1]
			RFF = header[3]
			ALT = header[4]
			SAMPLE = header[9:]
			SAMPLE = '\t'.join(SAMPLE)
			OUT1.write(CHR+'\t'+POS+'\t'+RFF+'\t'+ALT+'\t'+SAMPLE+'\n')
			OUT2.write(CHR+'\t'+POS+'\t'+RFF+'\t'+ALT+'\t'+SAMPLE+'\n')
			snp_count = defaultdict(int)
			snp_count2 = defaultdict(list)
			for i in f_csv:
				sample = ' '.join(i[9:])
				sample_tmp = re.findall(r'(./.)',sample)
				locate = i[3]+i[4]
				snp_count[locate] += 1 # count genotype

				snp_count2[locate].append(sample_tmp.count('0/0')) # count genotype in every sample
				snp_count2[locate].append(sample_tmp.count('0/1'))
				snp_count2[locate].append(sample_tmp.count('1/1'))
				########################################################
				geno = []
				if select == 'a': #simple
					for ii in sample_tmp:
						if ii == '0/0':
							geno.append(i[3])
						elif ii == '0/1':
							geno.append(pattern[locate])
						elif ii == '1/1':
							geno.append(i[4])
						elif ii == './.':
							geno.append('-')
						else:
							pass
				else: #double
					for ii in sample_tmp:
						if ii == '0/0':
							geno.append(i[3]+i[3])
						elif ii == '0/1':
							geno.append(i[3]+i[4])
						elif ii == '1/1':
							geno.append(i[4]+i[4])
						elif ii == './.':
							geno.append('--')
						else:
							pass
				sample_vcf = '\t'.join(sample_tmp)
				sample_geno = '\t'.join(geno)
				subject_vcf = i[0]+'\t'+i[1]+'\t'+i[3]+'\t'+i[4]+'\t'+sample_vcf+'\n'
				subject_geno = i[0]+'\t'+i[1]+'\t'+i[3]+'\t'+i[4]+'\t'+sample_geno+'\n'
				OUT1.write(subject_vcf)
				OUT2.write(subject_geno)
			#############################################################
			final_dict = {}
			for key,value in snp_count2.items():
				number_aa = 0
				number_Aa = 0
				number_AA = 0
				for count in range(0,len(value),3): # add aa Aa AA
					number_aa += value[count]
				for count in range(0,len(value),3):
					number_Aa += value[count+1]
				for count in range(0,len(value),3):
					number_AA += value[count+2]
				final_dict[key[0]+' --> '+key[1]]=[number_aa,number_Aa,number_AA]

			snp_type = snp_count.keys()
			snp_type = [i[0]+' --> '+i[1] for i in snp_type]
			snp_number = snp_count.values()
			snp_tuple = []
			for i in zip(snp_type,snp_number):
				snp_tuple.append(i)
			snp_tuple = sorted(snp_tuple, key=lambda x: x[1]) #sort by snp number
			snp_type = [i[0] for i in snp_tuple]
			snp_number = [i[1] for i in snp_tuple]
			heatmap_list = []
			for i in snp_type[::-1]:
				heatmap_list.append(final_dict[i])
			####################################### make pdf
			fig = plt.figure(figsize=(10,8))
			plt.style.use('ggplot')
			plt.subplots_adjust(wspace =0)
			ax = fig.add_subplot(1,2,2)
			plt.barh(range(len(snp_type)),snp_number,tick_label=list(snp_type),color=['royalblue'])
			plt.xlabel('snp number',{'family' : 'Arial','weight' : 'normal','size' : 12})
			ax.xaxis.get_major_formatter().set_powerlimits((0,1))
			ax.yaxis.set_ticks_position('right') 
			plt.tick_params(labelsize=12)

			plt.subplot(1,2,1)
			ax2 = fig.add_subplot(1,2,1)
			im = plt.imshow(heatmap_list, cmap=plt.cm.hot_r,interpolation='none')
			cbaxes = fig.add_axes([0.08, 0.75, 0.02, 0.2]) #set colorbar position
			xx = plt.colorbar(im,cbaxes)
			#plt.xticks(range(3), ['aa','Aa','AA'])
			plt.tick_params(labelsize=12)
			ax2.yaxis.set_ticks([])
			ax2.xaxis.set_ticks([])
			plt.tight_layout()
			plt.savefig("snp_count.pdf")
	except FileNotFoundError:
		print("The file you inputed was not found")
	except:
		print("Please input your option, or use --help")