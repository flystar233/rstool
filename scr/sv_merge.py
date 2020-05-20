import csv
import re
import itertools
import time
import click
def decomment(csvfile):
	for row in csvfile:
		raw = row.split('#')[0].strip()
		if raw: yield raw
def pos_filter(sv_data):
	##########################删除start位置是相同的行（只保留一个）
	sv_i12 = []
	for i in sv_data:
		sv_i12.append(i[0]+" "+str(i[1]))
	sv_i12_uniq = set(sv_i12)
	sv_i12_uniq_index = []
	for i in sv_i12_uniq:
		sv_i12_uniq_index.append(sv_i12.index(i))
	sv2= []
	for i in sorted(sv_i12_uniq_index):
		sv2.append(sv_data[i])
	#########################获得start位置pos之间距离小于1000bp点的索引
	sv2_1 = []
	sv2_2 = []
	for i in sv2:
		sv2_1.append(int(i[1]))
	for i in sv2:
		sv2_2.append(i[0])
	pos_500 = []
	for i in range(len(sv2_1)):
		if i >=len(sv2_1)-1:
			break
		else:
			cha = sv2_1[i+1]-sv2_1[i]
			if sv2_1[i] < 10000000 and cha>=-1000 and cha < 1000:
				pos_500.append(sv2_2[i]+" "+str(sv2_1[i]))
				pos_500.append(sv2_2[i+1]+" "+str(sv2_1[i+1]))
			elif sv2_1[i] > 10000000 and cha>=-2000 and cha < 2000:
				pos_500.append(sv2_2[i]+" "+str(sv2_1[i]))
				pos_500.append(sv2_2[i+1]+" "+str(sv2_1[i+1]))
			else:
				pass
	pos_500 =sorted(set(pos_500))
	##########################对小于500bp点过滤
	count =0
	final_sv =[]
	for i in sv2:
		if i[0]+" "+str(i[1]) in pos_500:
			pass
		else:
			final_sv.append(i)
			count+=1
	return final_sv
def merge_two(csv1,csv2,cut):
	tmp_list = []
	for bm in itertools.product(csv1,csv2):
		if bm[0][0]==bm[1][0] and bm[0][3]==bm[1][3]:
			head = max([int(bm[0][1]),int(bm[1][1])]) #获取大头
			tail = min([int(bm[0][2]),int(bm[1][2])]) #获取小尾
			cha = int(tail)-int(head) #获取overlap
			if  cha > 0 and cha/(int(bm[0][2])-int(bm[0][1])) >=cut and cha/(int(bm[1][2])-int(bm[1][1])) >=cut:
				tmp_list1 = (bm[0][0],bm[0][1],bm[0][2],bm[1][1],bm[1][2],cha/(int(bm[0][2])-int(bm[0][1]))+cha/(int(bm[1][2])-int(bm[1][1])),bm[0][5]+";"+bm[1][5])
				tmp_list.append(tmp_list1)
			else:
				pass
		else:
			pass
	return tmp_list
def uniq_sv(tmp_list,type_sv):
	tmp_list2 =sorted(tmp_list, key=lambda x:(x[0],int(x[1]),-x[5])) #-x[5] 对overlap所占2个片段百分比之和进行降序排序
	sv_i12 = []
	for i in tmp_list2:
		sv_i12.append(i[0]+" "+i[1])
	sv_i12_uniq = set(sv_i12) #取uniq
	sv_i12_uniq_index = []
	for i in sv_i12_uniq:
		sv_i12_uniq_index.append(sv_i12.index(i)) #index只会取匹配到的第一个列表元素
	sv_final_tmp= []
	for i in sorted(sv_i12_uniq_index):
		sv_final_tmp.append(tmp_list2[i])
	sv_final = []
	for i in sv_final_tmp:
		start = int((int(i[1])+int(i[3]))/2)
		end = int((int(i[2])+int(i[4]))/2)
		sv_ = (i[0],start,end,type_sv,end-start,i[6])
		sv_final.append(sv_)
	return sv_final
def merge_sv(method1,method2,cut,TYPE,out,method3=None):
	with open(method1,'r') as IN1,open(method2,'r') as IN2,open(out,'w',newline='') as OUT:
		header = ['#chromosome','start','end','type','size','info']
		b_csv = csv.reader(decomment(IN1),delimiter='\t')
		m_csv = csv.reader(decomment(IN2),delimiter='\t')
		if method3:
			with open(method3,'r') as IN3:
				s_csv = csv.reader(decomment(IN3),delimiter='\t')
				tmp_list = merge_two(b_csv,m_csv,cut) #使用接口
				sv_final = uniq_sv(tmp_list,type_sv=TYPE) #使用接口
				tmp_list_final = []
				for bm in itertools.product(sv_final,s_csv):
					if bm[0][0]==bm[1][0] and bm[0][3]==bm[1][3]:
						head = max([int(bm[0][1]),int(bm[1][1])])
						tail = min([int(bm[0][2]),int(bm[1][2])])
						cha = int(tail)-int(head)
						if  cha > 0 and cha/(int(bm[0][2])-int(bm[0][1])) >=cut and cha/(int(bm[1][2])-int(bm[1][1])) >=cut:
							start = int((int(bm[0][1])+int(bm[1][1]))/2)
							end = int((int(bm[0][2])+int(bm[1][2]))/2)
							sv_ = (bm[0][0],start,end,TYPE,end-start,bm[0][5]+";"+bm[1][5])
							tmp_list_final.append(sv_)
				bs = csv.writer(OUT,delimiter='\t')
				bs.writerow(header)
				bs.writerows(tmp_list_final)
		else:
			tmp_list = merge_two(b_csv,m_csv,cut)
			sv_final = uniq_sv(tmp_list,type_sv=TYPE)
			bs = csv.writer(OUT,delimiter='\t')
			bs.writerow(header)
			bs.writerows(sv_final)

@click.command("sv_merge", help="sv_merge by breakdancer,smoove and manta.")
@click.option("-b", "--breakdancer", type=click.Path(exists=True), help="breakdancer file")
@click.option("-m", "--manta", type=click.Path(exists=True), help="manta vcf file")
@click.option("-s", "--smoove", type=click.Path(exists=True),help="smoove vcf file")
@click.option(
    "-c", "--cutover", default=0.8, type=float, help="have reciprocal coordinate overlap larger than this threshold would be clustered together",show_default=True)
def sv_merge(breakdancer,manta,smoove,cutover):
	with open(breakdancer) as IN,open('breakdancer_DEL.txt','w',newline='') as b1,open('breakdancer_INV.txt','w',newline='') as b2:
		f_csv = csv.reader(decomment(IN),delimiter='\t') # skip annoation(#)
		breakdancer_header = ['#chromosome','start','end','type','size','info']
		row_breakdancer_DEL = []
		row_breakdancer_INV = []
		for sv in f_csv:
			tmp1 = []
			tmp2 = []
			if (abs(int(sv[7])) < 100 or int(sv[8]) < 60 or int(sv[9]) < 5):
				pass
			else:
				if sv[6]=='DEL':
					tmp1 = [sv[0],sv[1],sv[4],sv[6],sv[7],'score='+sv[8]+';'+'PE='+sv[9]]
					tmp_del =tuple(tmp1)
					row_breakdancer_DEL.append(tmp_del)
				elif sv[6]=='INV':
					tmp2 = [sv[0],sv[1],sv[4],sv[6],sv[7],'score='+sv[8]+';'+'PE='+sv[9]]
					tmp_inv =tuple(tmp2)
					row_breakdancer_INV.append(tmp_inv)
				else:
					pass
		final_DEL = pos_filter(row_breakdancer_DEL)
		final_INV = pos_filter(row_breakdancer_INV)
		b1_csv = csv.writer(b1,delimiter='\t')
		b1_csv.writerow(breakdancer_header)
		b1_csv.writerows(final_DEL)
		b2_csv = csv.writer(b2,delimiter='\t')
		b2_csv.writerow(breakdancer_header)
		b2_csv.writerows(final_INV)
		print(">>> End of breakdancer file processing.")
	with open(smoove) as IN,open('smoove_DEL.txt','w',newline='') as s1,open('smoove_INV.txt','w',newline='') as s2,open('smoove_DUP.txt','w',newline='') as s3:
		f_csv = csv.reader(decomment(IN),delimiter='\t') # skip annoation(#)
		smoove_header = ['#chromosome','start','end','type','size','info']
		row_smoove_DEL = []
		row_smoove_INV = []
		row_smoove_DUP = []
		for sv in f_csv:
			tmp1 = []
			tmp2 = []
			tmp3 = []
			if (re.findall(r'SVTYPE=BND;',sv[7]) or abs(int(re.findall(r'SVLEN=(-*\d+);',sv[7])[0])) < 100 or int(re.findall(r':(\d+)',sv[7])[0]) < 5):
				pass
			else:
				if re.findall(r'SVTYPE=DEL;',sv[7]):
					tmp1 = [sv[0],sv[1],re.findall(r'END=(\d+);',sv[7])[0],re.findall(r'SVTYPE=(\w+);',sv[7])[0],re.findall(r'SVLEN=(-*\d+);',sv[7])[0],'SU='+re.findall(r':(\d+)',sv[7])[0]]
					tmp_del =tuple(tmp1)
					row_smoove_DEL.append(tmp_del)
				elif re.findall(r'SVTYPE=INV;',sv[7]):
					tmp2 = [sv[0],sv[1],re.findall(r'END=(\d+);',sv[7])[0],re.findall(r'SVTYPE=(\w+);',sv[7])[0],re.findall(r'SVLEN=(-*\d+);',sv[7])[0],'SU='+re.findall(r':(\d+)',sv[7])[0]]
					tmp_inv =tuple(tmp2)
					row_smoove_INV.append(tmp_inv)
				elif re.findall(r'SVTYPE=DUP;',sv[7]):
					tmp3 = [sv[0],sv[1],re.findall(r'END=(\d+);',sv[7])[0],re.findall(r'SVTYPE=(\w+);',sv[7])[0],re.findall(r'SVLEN=(-*\d+);',sv[7])[0],'SU='+re.findall(r':(\d+)',sv[7])[0]] 
					tmp_dup =tuple(tmp3)
					row_smoove_DUP.append(tmp_dup)
				else:
					pass
		final_DEL = pos_filter(row_smoove_DEL)
		final_INV = pos_filter(row_smoove_INV)
		final_DUP = pos_filter(row_smoove_DUP)
		s1_csv = csv.writer(s1,delimiter='\t')
		s1_csv.writerow(smoove_header)
		s1_csv.writerows(final_DEL)
		s2_csv = csv.writer(s2,delimiter='\t')
		s2_csv.writerow(smoove_header)
		s2_csv.writerows(final_INV)
		s3_csv = csv.writer(s3,delimiter='\t')
		s3_csv.writerow(smoove_header)
		s3_csv.writerows(final_DUP)
		print(">>> End of smoove file processing.")
	with open(manta) as IN,open('manta_DEL.txt','w',newline='') as m1,open('manta_INV.txt','w',newline='') as m2,open('manta_DUP.txt','w',newline='') as m3:
		f_csv = csv.reader(decomment(IN),delimiter='\t') # skip annoation(#)
		manta_header = ['#chromosome','start','end','type','size','info']
		row_manta_DEL = []
		row_manta_INV = []
		row_manta_DUP = []
		for sv in f_csv:
			tmp1 = []
			tmp2 = []
			tmp3 = []
			if ( sv[6]== 'PASS' and (sv[4]=='<DEL>' or sv[4]=='<DUP:TANDEM>' or sv[4]=='<INV>')):
				if re.findall(r'SVTYPE=DEL',sv[7]) and abs(int(re.findall(r'SVLEN=(-*\d+)',sv[7])[0])) > 100:
					tmp1 = [sv[0],sv[1],re.findall(r'END=(\d+)',sv[7])[0],re.findall(r'SVTYPE=(\w+)',sv[7])[0],re.findall(r'SVLEN=(-*\d+)',sv[7])[0],'.']
					tmp_del =tuple(tmp1)
					row_manta_DEL.append(tmp_del)
				elif re.findall(r'SVTYPE=INV',sv[7]) and abs(int(re.findall(r'SVLEN=(-*\d+)',sv[7])[0])) > 100:
					tmp2 = [sv[0],sv[1],re.findall(r'END=(\d+)',sv[7])[0],re.findall(r'SVTYPE=(\w+)',sv[7])[0],re.findall(r'SVLEN=(-*\d+)',sv[7])[0],'.']
					tmp_inv =tuple(tmp2)
					row_manta_INV.append(tmp_inv)
				elif re.findall(r'SVTYPE=DUP',sv[7]) and abs(int(re.findall(r'SVLEN=(-*\d+)',sv[7])[0])) > 100:
					tmp3 = [sv[0],sv[1],re.findall(r'END=(\d+)',sv[7])[0],re.findall(r'SVTYPE=(\w+)',sv[7])[0],re.findall(r'SVLEN=(-*\d+)',sv[7])[0],'.']
					tmp_dup =tuple(tmp3)
					row_manta_DUP.append(tmp_dup)
			else:
				pass
		final_DEL = pos_filter(row_manta_DEL)
		final_INV = pos_filter(row_manta_INV)
		final_DUP = pos_filter(row_manta_DUP)
		m1_csv = csv.writer(m1,delimiter='\t')
		m1_csv.writerow(manta_header)
		m1_csv.writerows(final_DEL)
		m2_csv = csv.writer(m2,delimiter='\t')
		m2_csv.writerow(manta_header)
		m2_csv.writerows(final_INV)
		m3_csv = csv.writer(m3,delimiter='\t')
		m3_csv.writerow(manta_header)
		m3_csv.writerows(final_DUP)
		print(">>> End of manta file processing.")
	print(">>> Start merging files of the same variant type")
	merge_sv('breakdancer_DEL.txt','manta_DEL.txt',cutover,'DEL','breakdancer_manta_smoove_DEL.txt','smoove_DEL.txt')
	merge_sv('breakdancer_INV.txt','smoove_INV.txt',cutover,'INV','breakdancer_smoove_INV.txt')
	merge_sv('smoove_DUP.txt','manta_DUP.txt',cutover,'DUP','smoove_manta_DUP.txt')
	print(">>> All done!")
