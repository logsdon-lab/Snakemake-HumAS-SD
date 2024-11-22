# Convert monomeric bed file into the HOR bed file
# Usage: python3 hor2stv.py <ASat live HORs>.bed > <StV>.bed
import argparse


def stv_namer(live_stv_name, mons_numbers, strand):
    if strand == '+':
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        else: # 8&12 in chr18
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) + 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            # 8&12 in chr18
            else:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    # reversed. strand == '-'
    else:
        i_prev = 25
        stv_name = mons_numbers[0]
        if mons_numbers[0].isdigit():
            i_prev = mons_numbers[0]
        elif '/' in mons_numbers[0] and 'S' not in mons_numbers[0]: # first mon is hybrid
            i_prev = 25
        status = 'Closed'
        for i in mons_numbers[1:]:
            if i.isdigit():
                if int(i) == int(i_prev) - 1:
                    i_prev = i
                    status = 'ToBeClosed'
                else:
                    if status == 'ToBeClosed':
                        stv_name += '-{}_{}'.format(i_prev, i)
                    else:
                        stv_name += '_{}'.format(i)
                    i_prev = i
                    status = 'Closed'
            # hybrid like '4/7'
            elif '/' in i and 'S' not in i:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                i_prev = 25
                status = 'Closed'
            else:
                if status == 'ToBeClosed':
                    stv_name += '-{}_{}'.format(i_prev, i)
                else:
                    stv_name += '_{}'.format(i)
                status = 'Closed'
                i_prev = '20'
        if status == 'ToBeClosed':
            stv_name += '-{}'.format(i)
    stv_name = '{}.{}'.format(live_stv_name, stv_name)
    return(stv_name)

def print_bed(bed,outfile):
    outfile = open(outfile,'w')
    for i in bed:
        # [contig, hor_start, hor_end, stv_name, '0', strand, hor_start, hor_end, '0,0,0']
        outfile.write(i[0]+'\t' + i[1] + '\t' + i[2] + '\t' + i[3] + '\t' + i[4] + '\t' + i[5] + '\t' + i[6] + '\t' + i[7]  + '\t' + i[8]+'\n')
    outfile.close()

def main():

    parser = argparse.ArgumentParser(description="Annotation HOR for centromere reads")
    parser.add_argument("-i", "--input_bed_path", help="", required=True)
    parser.add_argument("-o", "--output_bed_path", help="", required=True)
    args = parser.parse_args()

    input_bed_path = args.input_bed_path
    output_bed_path = args.output_bed_path

    # NA19650_rc-chr22_h1tg000022l#1-28700957:2895006-6578064_renamed.bed
    # input_bed_path = 'D:/working/HOR_STV/test_case/input_horname_test.bed'
    # output_bed_path = 'D:/working/HOR_STV/test_case/output_horname_test.bed'

    # parse input bed, #remove NOT HOR
    input_bed = []
    with open(input_bed_path) as bed:
        for line in bed:
            if line[:5] != 'track': # skip header
                monomer_name = line.split('\t')[3]
                contig_name = line.split('\t')[0]
                if contig_name == 'chr1' or contig_name == 'chr19':
                    if 'S1C1/5/19H1L.6/4' in line:
                        line = line.replace('S1C1/5/19H1L.6/4', 'S1C1/5/19H1L.6')

                if contig_name == 'chr5':
                    if 'S1C1/5/19H1L.2/6' in line:
                        line = line.replace('S1C1/5/19H1L.2/6', 'S1C1/5/19H1L.6')

                if contig_name == 'chr8':
                    if 'S2C8H1L.6/7s' in line:
                        line = line.replace('S2C8H1L.6/7s', 'S2C8H1L.6/7')

                if contig_name == 'chr13':
                    if 'S2C13/21H1-B.10' in line:
                        line = line.replace('S2C13/21H1-B.10', 'S2C13/21H1L.10')

                if contig_name == 'chr18':
                    if 'S2C18H2-E.x/2'  in line:
                        line = line.replace('S2C18H2-E.x/2', 'S2C18H2-E.2')
                    if 'S2C18H2-E.2/x'  in line:
                        line = line.replace('S2C18H2-E.2/x', 'S2C18H2-E.2')

                if monomer_name.startswith('S'):
                    input_bed.append(line.split())

    # list of contigs
    contigs = []
    for line in input_bed:
        if line[0] not in contigs:
            contigs.append(line[0])

    # MAIN LOOP
    stv_bed = []
    for contig in contigs:
        #print(contig)
        live_mons = []
        for line in input_bed:
            if line[0] == contig:
                live_mons.append(line)

        # 在monomer seq中检测断裂点

        # splite monomer_seq based on different +-, gap and different HOR name
        monomer_seqs = []
        seq = [live_mons[0]]
        for index in range(len(live_mons) - 1):
            next = index + 1
            last = index
            next_hor_name = live_mons[next][3].split('.')[0]
            last_hor_name = live_mons[last][3].split('.')[0]

            if int(live_mons[next][1]) - int(seq[-1][2]) > 160: # gaps
                monomer_seqs.append(seq)
                seq = [live_mons[next]]
            elif live_mons[next][5] != seq[-1][5]: # different +-
                monomer_seqs.append(seq)
                seq = [live_mons[next]]
            elif next_hor_name != last_hor_name:
                monomer_seqs.append(seq)
                seq = [live_mons[next]]
            else:
                seq.append(live_mons[next])
        monomer_seqs.append(seq)
        HOR_units = []
        for seq in monomer_seqs:
            # detect break points
            if len(seq) == 1:
                HOR_units.append(seq)
                continue
            start_index = 0
            first_monomer = seq[0]
            strand = first_monomer[5]
            if strand == '+':
                # break points is  5 4, last >= next 5 5
                for m in range(len(seq) - 1):
                    last = seq[m]
                    next = seq[m+1]
                    monomer_last = last[3].split('.')[-1]
                    monomer_next = next[3].split('.')[-1]
                    monomer_id_last = ''
                    monomer_id_next = ''
                    if '/' in monomer_last:
                        monomer_id_last = int(monomer_last.split('/')[-1])
                    else:
                        monomer_id_last = int(monomer_last)
                    if '/' in monomer_next:
                        monomer_id_next = int(monomer_next.split('/')[0])
                    else:
                        monomer_id_next = int(monomer_next)

                    if monomer_id_next < monomer_id_last:
                        HOR_units.append(seq[start_index:m + 1])
                        start_index = m + 1


            else:
                # break points is 5 4 3 * 5 4 3 2 1, last <= next 5 5
                for m in range(len(seq) - 1):
                    last = seq[m]
                    next = seq[m+1]
                    monomer_last = last[3].split('.')[-1]
                    monomer_next = next[3].split('.')[-1]
                    monomer_id_last = ''
                    monomer_id_next = ''
                    if '/' in monomer_last:
                        monomer_id_last = int(monomer_last.split('/')[-1])
                    else:
                        monomer_id_last = int(monomer_last)
                    if '/' in monomer_next:
                        monomer_id_next = int(monomer_next.split('/')[0])
                    else:
                        monomer_id_next = int(monomer_next)

                    if monomer_id_next > monomer_id_last:
                        HOR_units.append(seq[start_index:m + 1])
                        start_index = m + 1

            if len(seq[start_index:len(seq)]) != 0:
                HOR_units.append(seq[start_index:len(seq)])

        stvs = []
        for hor in HOR_units:
            hor_strand = hor[0][5]
            name = hor[0][3].split('.')[0]
            hor_start = hor[0][1]
            hor_end = hor[-1][2]
            monomer_seq = []
            for m in hor:
                monomer_seq.append(m[3].split('.')[-1])
            stv_name = stv_namer(name, monomer_seq, hor_strand)
            stvs.append([contig, hor_start, hor_end, stv_name, '0', hor_strand, hor_start, hor_end, '0,0,0'])

        stv_bed += stvs

    print_bed(stv_bed,output_bed_path)


if __name__ == '__main__':
    main()
