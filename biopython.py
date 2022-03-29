from datetime import datetime
from xml.dom.pulldom import PROCESSING_INSTRUCTION
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio import SearchIO
import sys
import pandas as pd
import os
from jinja2 import Environment, FileSystemLoader
import pdfkit
import copy
# cmd:
# python biopython.py db/Omicron.fasta 0
if len(sys.argv) != 3:
    print('usage : python biopython.py db/Omicron.fasta 0 ')
    exit(0)

if os.path.exists('origin_data/') == False:
    os.mkdir('origin_data')
if os.path.exists('result/') == False:
    os.mkdir('result')
if os.path.exists('db/') == False:
    os.mkdir('db')

db_filename = sys.argv[1]
# 建库
if int(sys.argv[2]) == 100:
    local_db = NcbimakeblastdbCommandline(
        dbtype="nucl", input_file=db_filename)()[0]
    exit(0)


type_index = int(sys.argv[2])
filename_prefixs = ['N_', 'O_']
fname = os.path.splitext(os.path.basename(db_filename))[0]
filename = filename_prefixs[type_index] + fname
output_xml_filename = 'origin_data/' + filename+'.xml'
report_filename = 'result/'+filename + '_report.html'
detail_html = 'origin_data/' + filename+'.html'
detail_pdf = 'result/' + filename+'.pdf'
typenames = ['N', 'ORF1ab']

df = pd.read_excel('database.xlsx', sheet_name=0)
database = {}
for line in df.index:
    temp = {}
    temp['duzhuname'] = df[df.columns[1]].loc[line]
    pair = {
        'B_1_427_B_1_429': 'B.1.427_B.1.429'
    }
    if temp['duzhuname'] in pair:
        temp['duzhuname'] = pair[temp['duzhuname']]

    typename = df[df.columns[2]].loc[line].strip()
    if typenames[type_index] != typename:
        continue
    temp['valid_id'] = df[df.columns[5]].loc[line]
    temp['valid_time'] = str(df[df.columns[4]].loc[line])
    if '00:00' in temp['valid_time']:
        temp['valid_time'] = temp['valid_time'].split(' ')[0]

    temp['description'] = str(df[df.columns[6]].loc[line]).strip()
    temp['yingxiang'] = str(df[df.columns[9]].loc[line])
    temp['yanzhengjieguo'] = str(df[df.columns[10]].loc[line])
    temp['suoyouba'] = str(df[df.columns[11]].loc[line])

    if temp['yingxiang'].upper() == 'NAN':
        temp['yingxiang'] = ''
    if temp['yanzhengjieguo'] .upper() == 'NAN':
        temp['yanzhengjieguo'] = ''
    if temp['suoyouba'] .upper() == 'NAN':
        temp['suoyouba'] = ''
    if temp['description'] in database:
        if fname in temp['duzhuname']:
            database[temp['description']] = temp
        continue
    database[temp['description']] = temp

query_filename = ['N_query.fasta', 'O_query.fasta']
query = SeqIO.read(query_filename[type_index], 'fasta')

db_size = len(list(SeqIO.parse(db_filename, 'fasta')))
cline = NcbiblastnCommandline(query=query_filename[type_index],
                              line_length=len(query.seq),
                              max_target_seqs=db_size,
                              db=db_filename,
                              evalue=1,
                              out=output_xml_filename,
                              html=True, outfmt=5)()

blast_qresults = SearchIO.read(output_xml_filename, "blast-xml")
if len(blast_qresults) != db_size:
    print('len(blast_qresults) != db_size', len(blast_qresults), db_size)


ans = [
    'TCAACTCCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTGGTAAAGGCCAACAACAACAA',
    'CCTACAACTTGTGCTAATGACCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTG'
]

ans_html = [
    'TCAACTCCAGGCAGCAGTA<span style="background:#90ee90">GGGGAACTTCTCCTGCTAGAAT</span>GGCTGGCAATGGCGGTGATGCTGCTCTTGCT<span style="background:#90ee90">TTGCTGCTGCTTGACAGATT</span>GAAC<span style="background:#90ee90">CAGCTTGAGAGCAAAATGTCTG</span>GTAAAGGCCAACAACAACAA',
    'CCTACAACTTGTGCTAATGA<span style="background:#90ee90">CCCTGTGGGTTTTACACTTAA</span>AAACACAGTCTGTA<span style="background:#90ee90">CCGTCTGCGGTATGTGGAAAGGTTATGG</span>CTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAG<span style="background:#90ee90">TCAGCTGATGCACAATCGT</span>TTTTAAACGGGTTTGCGGTG'
]

positions = [
    # 19,22,31,20,4,22,20
    [0, 19, 41, 72, 92, 96, 118, 138],
    # 20,21,14,28,37,19,20
    [0, 20, 41, 55, 83, 120, 139, 159]
]

# 处理有插入情况下，引物探针的index


def get_key_pos_index(insert_indexs):
    pos = copy.deepcopy(positions[type_index])
    for index in insert_indexs:
        print(index)
        start = index[0]
        length = index[1] - index[0]
        for i in range(0, len(pos)):
            if pos[i] > start:
                pos[i] = pos[i] + length
    print(pos)
    return pos


def change_to_diff_html(db, query):
    retstr = ""
    if len(db) != len(query):
        print(db, query)
    if db == query:
        return '<span style="background:#90ee90">{}</span>'.format(query)
    for i in range(0, len(db)):
        # deal with insert
        if db[i] == '-':
            retstr += '<span style="background:#00ffff">{}</span>'.format(
                query[i])
        elif db[i] != query[i]:
            retstr += '<span style="background:red">{}</span>'.format(query[i])
        else:
            retstr += '<span style="background:#90ee90">{}</span>'.format(
                query[i])
    return retstr


def get_key_compstr(seq, pos):
    return seq[pos[1]:pos[2]] + seq[pos[3]:pos[4]]+seq[pos[5]:pos[6]]


def change_all_seq_to_html(query, seq, pos):
    #pos = positions[type_index]
    s1, s2, s3, s4 = seq[pos[0]:pos[1]],  seq[pos[2]
        :pos[3]], seq[pos[4]:pos[5]], seq[pos[6]:]

    # 正向引物,探针，反向引物
    ss1, ss2, ss3 = seq[pos[1]:pos[2]], seq[pos[3]:pos[4]], seq[pos[5]:pos[6]]

    r1 = change_to_diff_html(query[pos[1]:pos[2]], seq[pos[1]:pos[2]])
    r2 = change_to_diff_html(query[pos[3]:pos[4]], seq[pos[3]:pos[4]])
    r3 = change_to_diff_html(query[pos[5]:pos[6]], seq[pos[5]:pos[6]])
    return s1+r1+s2+r2+s3+r3+s4


def get_diff_range(hit, query):
    ret = []
    if len(hit) != len(query):
        print('str len error')
        exit(0)
    i = 0
    insert_count = 0
    while i < len(hit):
        if hit[i] != query[i]:
            temp = []
            is_deletion = False
            is_insertion = False
            if hit[i] == '-':
                is_deletion = True
            elif query[i] == '-':
                is_insertion = True
                insert_count = insert_count
            temp.append(i - insert_count)
            j = i+1
            while j < len(hit):
                if hit[j] == query[j]:
                    break
                if is_deletion and hit[j] != '-':
                    break
                if is_deletion == False and hit[j] == '-':
                    break
                if is_insertion and query[i] != '-':
                    break
                if is_insertion == False and query[i] == '-':
                    break
                j = j + 1
            temp.append(j-insert_count)
            if is_deletion:
                temp.append('-')
            elif is_insertion:
                temp.append('+')
            ret.append(temp)
            i = j - 1
        i = i+1
    return ret


def des_helper(t1, t2, typestr, reserve):
    if t1 == t2:
        return ""
    if reserve:
        t1 = t1[::-1]
        t2 = t2[::-1]
    retstr = typestr
    diffs = get_diff_range(t1, t2)

    for diff in diffs:
        if diff[1] - diff[0] == 1:
            if len(diff) == 3 and diff[2] == '-':
                retstr += '3\'端第{}个碱基删除{}突变'.format(
                    diff[0]+1, t2[diff[0]:diff[1]])
            elif len(diff) == 3 and diff[2] == '+':
                retstr += '3\'端第{}个碱基插入{}突变'.format(
                    diff[0], t1[diff[0]:diff[1]])
            else:
                retstr += '3\'端第{}个碱基{}>{}突变'.format(
                    diff[0]+1, t2[diff[0]:diff[1]], t1[diff[0]:diff[1]])
        else:
            if reserve:
                if len(diff) == 3 and diff[2] == '-':
                    retstr += '3\'端第{}_{}个碱基删除{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]][::-1])
                elif len(diff) == 3 and diff[2] == '+':
                    retstr += '3\'端第{}个碱基插入{}突变'.format(
                        diff[0], t1[diff[0]:diff[1]][::-1])
                else:
                    retstr += '3\'端第{}_{}个碱基{}>{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]][::-1], t1[diff[0]:diff[1]][::-1])
            else:
                if len(diff) == 3 and diff[2] == '-':
                    retstr += '3\'端第{}_{}个碱基删除{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]])
                elif len(diff) == 3 and diff[2] == '+':
                    retstr += '3\'端第{}个碱基插入{}突变'.format(
                        diff[0], t1[diff[0]:diff[1]])
                else:
                    retstr += '3\'端第{}_{}个碱基{}>{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]], t1[diff[0]:diff[1]])
    retstr = retstr + ','
    return retstr


def des_helper2(t1, t2, typestr, reserve):
    if t1 == t2:
        return ""
    if reserve:
        t1 = t1[::-1]
        t2 = t2[::-1]
    retstr = typestr
    diffs = get_diff_range(t1, t2)
    for diff in diffs:
        if diff[1] - diff[0] == 1:
            if len(diff) == 3 and diff[2] == '+':
                retstr += '3\'端第{}个碱基插入{}突变'.format(
                    diff[0], t1[diff[0]:diff[1]][::-1])
            else:
                retstr += '3\'端第{}个碱基{}>{}突变'.format(
                    diff[0]+1, t2[diff[0]:diff[1]], t1[diff[0]:diff[1]])
        else:
            if reserve:
                if len(diff) == 3 and diff[2] == '+':
                    retstr += '3\'端第{}个碱基插入{}突变'.format(
                        diff[0], t1[diff[0]:diff[1]][::-1])
                else:
                    retstr += '3\'端第{}_{}个碱基{}>{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]][::-1], t1[diff[0]:diff[1]][::-1])
            else:
                if len(diff) == 3 and diff[2] == '+':
                    retstr += '3\'端第{}个碱基插入{}突变'.format(
                        diff[0], t1[diff[0]:diff[1]])
                else:
                    retstr += '3\'端第{}_{}个碱基{}>{}突变'.format(
                        diff[0]+1, diff[1], t2[diff[0]:diff[1]], t1[diff[0]:diff[1]])
    retstr = retstr + ','
    return retstr


def generate_desc_str(seq, answer, key_position):
    l1 = list(seq)

    pos = key_position
    # 正向引物,探针，反向引物
    ss1, ss2, ss3 = seq[pos[1]:pos[2]], seq[pos[3]:pos[4]], seq[pos[5]:pos[6]]
    ans1, ans2, ans3 = answer[pos[1]:pos[2]
                              ], answer[pos[3]:pos[4]], answer[pos[5]:pos[6]]

    retstr = ""

    retstr += des_helper(ss1, ans1, '正向引物', True)
    retstr += des_helper(ss3, ans3, '反向引物', False)
    retstr += des_helper(ss2, ans2, '探针', True)
    if len(retstr) == 0:
        return "", ""
    if retstr[-1] == ',':
        retstr = retstr[:-1]
    retstr2 = retstr
    retstr = ""
    retstr += des_helper2(ss1, ans1, '正向引物', True)
    retstr += des_helper2(ss3, ans3, '反向引物', False)
    retstr += des_helper2(ss2, ans2, '探针', True)
    if len(retstr) == 0:
        return ""
    if retstr[-1] == ',':
        retstr = retstr[:-1]

    return retstr, retstr2


def need_experiment(s, answer, pos):
    keystr = get_key_compstr(s, pos)
    strans = get_key_compstr(answer, pos)
    # 缺失多余5个，有N都不验证
    if keystr.count('-') >= 5:
        return False
    if 'N' in keystr:
        return False
    # 存在碱基的插入、缺失
    if '-' in answer or '-' in keystr:
        return True
    # 探针/引物上的突变加起来超过2个
    diff_count = 0
    for i in range(0, len(strans)):
        if strans[i] != keystr[i]:
            diff_count = diff_count + 1
    if diff_count >= 2:
        return True

    # 探针/引物3’端5个碱基内有突变
    check_indexs = [[pos[2]-5, pos[2]], [pos[3]-5, pos[3]], [pos[5], pos[5]+5]]
    for check_index in check_indexs:
        for i in range(check_index[0], check_index[1]):
            if answer[i] != s[i]:
                return True
    return False
    # 中国上传
    # 突变频率> 0.1 %


def process_insert(str1, str2):
    ret = []
    i = 0
    while i < len(str1):
        if str1[i] == '-':
            temp = []
            temp.append(i)
            j = i+1
            while j < len(str1):
                if str1[j] != '-':
                    break
                j = j + 1
            temp.append(j)
            temp.append(str2[i:j])
            ret.append(temp)
            i = j
        i = i+1
    processed_str = ""
    for i in range(0, len(str1)):
        if str1[i] != '-':
            processed_str = processed_str+str2[i]
    #return ret, processed_str
    return ret


def generate_detail_html(report):

    env = Environment(loader=FileSystemLoader('./'))
    template = env.get_template('detail_template.html')
    with open(detail_html, 'w+', encoding='utf-8') as fout:
        html_content = template.render(
            ans_html=ans_html[type_index],
            items=report)
        fout.write(html_content)


def change_html_to_pdf(htmlfile, pdffile):
    wkhtmltopdf_path = './bin/wkhtmltopdf.exe'
    config = pdfkit.configuration(wkhtmltopdf=wkhtmltopdf_path,
                                  )
    # 生成pdf文件，to_file为文件路径
    pdfkit.from_file(htmlfile, pdffile, configuration=config)


lists = []
need_lists = []
need_done_count = 0
unneed_lists = []


visited = {}
report = {}
pos = positions[type_index]
items = []
for blast_result in blast_qresults:
    fragment = blast_result[0][0]
    from_china = False

    keys = fragment.hit.description.strip().split('|')
    if len(keys) == 3:
        key = fragment.hit.description.strip().split('|')[1]
        from_china = 'CHINA' in fragment.hit.description.strip().upper()
    else:
        key = fragment.hit.id.strip().split('|')[1]
        from_china = 'CHINA' in fragment.hit.id.strip().upper()
    key = key.strip()
    if key in visited:
        continue
    visited[key] = True

    if 'EPI_ISL' not in key:
        print('warning : EPI_ISL not in key')
    origin_seq = fragment.hit.seq.encode("ascii").decode("utf-8")
    seq = origin_seq
    origin_answer = fragment.query.seq.encode("ascii").decode("utf-8")

    insert_indexs = []
    key_position = pos
    if len(origin_answer) > len(ans[type_index]):
        insert_indexs = process_insert(origin_answer, seq)
        # print(insert_indexs,seq)
    if len(origin_answer) < len(ans[type_index]):
        # 如果有缺失，则补全
        origin_answer = ans[type_index]
        for i in range(0, fragment.query_range[0]):
            seq = '-' + seq
        for i in range(fragment.query_range[1], len(origin_answer)):
            seq = seq + '-'
        origin_seq = seq

    if len(insert_indexs) > 0:
        key_position = get_key_pos_index(insert_indexs)
    value = get_key_compstr(origin_seq, key_position)
    # 为了保证html和pdf中对齐，补空格
    for _ in range(len(key), 18):
        key = key + ' '

    oneitem = {}
    oneitem['key'] = key
    oneitem['start_index'] = fragment.hit_range[0]+1
    oneitem['htmlseq'] = change_all_seq_to_html(
        origin_answer, origin_seq, key_position)
    items.append(oneitem)
    if value in report:
        if from_china:
            report[value]['from_china'] = '是'
            report[value]['key'] = key
            report[value]['start_index'] = fragment.hit_range[0]+1
        report[value]['count'] = report[value]['count'] + 1
        continue
    desc_str, desc_str2 = generate_desc_str(
        origin_seq, origin_answer, key_position)
    report[value] = {}
    if from_china:
        report[value]['from_china'] = '是'
    else:
        report[value]['from_china'] = '否'
    report[value]['desc_str'] = desc_str
    report[value]['desc_str2'] = desc_str2
    report[value]['description'] = desc_str2
    report[value]['count'] = 1
    report[value]['key'] = key
    report[value]['seq'] = origin_seq
    report[value]['htmlseq'] = oneitem['htmlseq']
    report[value]['start_index'] = fragment.hit_range[0]+1
    report[value]['key_position'] = copy.deepcopy(key_position)
    report[value]['need_experiment'] = need_experiment(
        origin_seq, origin_answer, key_position)
    if from_china:
        report[value]['need_experiment'] = True

generate_detail_html(items)
change_html_to_pdf(detail_html, detail_pdf)

for key in report:
    # 去掉没有发生突变的结果
    if report[key]['description'] != "":
        lists.append(report[key])
lists = sorted(lists, key=lambda x: (-x['count'], x['key']))

for item in lists:
    ratio_str = '{:4.2f}%'.format(100.0 * item['count'] / len(items))
    item['ratio'] = ratio_str
    matched = False
    for t in database:
        if item['desc_str'] != database[t]['description'] and item['description'] != database[t]['description']:
            continue
        matched = True
        #item['description']  = desc_str2
        item['valid_id'] = database[t]['valid_id']
        item['valid_time'] = database[t]['valid_time']
        item['yingxiang'] = database[t]['yingxiang']
        item['yanzhengjieguo'] = database[t]['yanzhengjieguo']
        item['suoyouba'] = database[t]['suoyouba']

        item['status'] = '之前已验证'
        need_done_count = need_done_count+1
        need_lists.append(item)
        break
    if matched == False:
        if item['need_experiment'] == False:
            item['status'] = '无须验证'
            unneed_lists.append(item)
        else:
            item['status'] = '是'
            need_lists.append(item)


def generate_html(need, unneed, need_done_count):
    env = Environment(loader=FileSystemLoader('./'))
    template = env.get_template('report_template.html')
    with open(report_filename, 'w+', encoding='utf-8') as fout:
        html_content = template.render(
            need_done_count=need_done_count,
            total_count=len(items),
            ans_html=ans_html[type_index],
            need=need,
            unneed=unneed)
        fout.write(html_content)


generate_html(need_lists, unneed_lists, need_done_count)

result_filename = 'result/result.csv'
if os.path.exists(result_filename) == False:
    with open(result_filename, 'a+', encoding='gb18030') as fout:
        fout.write('突变类型,总突变序列数,总突变类型,需验证突变株,已验证,需要验证\n')
with open(result_filename, 'a+') as fout:
    fout.write('{},{},{},{},{},{}\n'.format(filename, len(items), len(unneed_lists) + len(need_lists),
                                            len(need_lists), need_done_count, len(need_lists) - need_done_count))

result_filename = 'result/detail.csv'
if os.path.exists(result_filename) == False:
    with open(result_filename, 'a+', encoding='gb18030') as fout:
        fout.write('毒株,基因,描述,数量,代表序列ID,是否需要验证,验证编号,验证时间,中国序列,碱基序列, \
        是否有影响,验证结果,突变株是否对所有靶标产生影响\n')
with open(result_filename, 'a+') as fout:
    for item in need_lists:
        if item['seq'][0] == '-':
            item['seq'] = ' ' + item['seq']
        fout.write('{},{},\"{}\",{},{},{},{},{},{},\"{}\",{},{},{}\n'.format(fname, typenames[type_index], item['description'],
                                                                             item['count'], item['key'].strip(), item['status'], item.get(
                                                                                 'valid_id', ""), item.get('valid_time', ""), item['from_china'], item['seq'],
                                                                             item.get('yingxiang', ""), item.get('yanzhengjieguo', ""), item.get('suoyouba', "")))
