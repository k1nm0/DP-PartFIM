#!/usr/bin/python3
# author liuxinyu
import math
import time
from collections import Counter
import numpy as np

def data_read(infile, t_c):
    f = open(infile, 'r', encoding='utf8')
    data = []  # [data_i,data_i,...]
    data_i = {}  # {item:{tran1,trans2,...}}
    trans = 0  # 事务数量
    trans_i = 0  # 分区事务数量
    times = 0  # 分区数
    maxtrans = []  # 存储分区中的事务长度
    MaxTrans = []  # 每个分区中位事务长度
    for row in f:
        # l=row.rstrip("\n").split(' ') #以“\n”结尾
        # l = row.rstrip(" \n").split(' ')# 以“ \n”结尾
        l=row.replace('\n','').replace(' ','')
        l=l.split('\t')
        trans += 1
        trans_i += 1
        for i in l:
            if i not in data_i:
                data_i[i] = set()
            data_i[i].add(trans)
            maxtrans.append(len(l))
        if trans_i > t_c:
            temp = data_i.copy()
            data.append(temp)
            trans_i = 0
            times += 1
            data_i.clear()
            L = sorted(maxtrans)
            temp_m = [L[len(L) // 2] if len(L) % 2 == 1 else "%.1f" % (0.5 * (L[len(L) // 2 - 1] + L[len(L) // 2])),
                      L[-1]]
            MaxTrans.append(temp_m)  # 将分区中位长度和最大长度存入MaxTrans
            maxtrans.clear()
    if trans_i > 0:
        data.append(data_i)
        L = sorted(maxtrans)
        temp_m = [L[len(L) // 2] if len(L) % 2 == 1 else "%.1f" % (0.5 * (L[len(L) // 2 - 1] + L[len(L) // 2])), L[-1]]
        MaxTrans.append(temp_m)
    return data, trans, times, trans_i, MaxTrans


def empty_file(outfile):
    with open(outfile, 'w') as f_output:
        f_output.close()
    file = open(outfile, "r+")
    file.truncate(0)
    file.close()


def write_into_file(string, outfile):
    with open(outfile, 'a') as f_output:
        f_output.write(string + "\n")
        f_output.close()



# 项集及其对应交易项和支持度
class Item():
    def __inti__(self, item, transet, sup):
        self.item = item
        self.transet = transet
        self.sup = sup


# 分区数据和分区支持度组合
class SupToData():
    def __init__(self, partsup, partdata, part):
        self.minsup = partsup
        self.data = partdata
        self.part = part


# 项集存在的数据分区
class item_part():
    def __init__(self, item):
        self.item = item
        self.data_i = []

    def count(self, data_i):
        self.data_i.append(data_i)


def data_sup(i):
    node = SupToData(partsup=part_minsup[i], partdata=data[i], part=i)
    return node


def partMine(sup_data):
    f_1 = []  # 元素：item,tranid,sup，频繁一项集
    F = {}  # item:item(item,transet,sup) ，存储最终的频繁项集
    Part_minsup = sup_data.minsup
    Part_data = sup_data.data
    Part_loc = sup_data.part
    # i_itemset(Part_data,Part_minsup,f_1,F)
    item_t = list(Part_data.items())
    for i in item_t:
        one_itemset(Part_data, Part_minsup, f_1, F, i, Part_loc)
    flist = sorted(f_1, key=lambda x: x[2], reverse=False)
    flist = [x[0] for x in flist]  # 频繁一项集item
    k_Mine({'0': set(flist)}, flist, F, Part_minsup, Part_loc)
    return {eval(x): F[x].sup for x in F.keys()}



def one_itemset(data, part_minsup, f_1, F, item_t, part):
    (item, t) = item_t
    i_sup = len(t)
    if i_sup < part_minsup:
        del data[item]
    else:
        f_1.append([item, t, i_sup])
        temp = Item()
        temp.item = item
        temp.transet = t
        temp.sup = i_sup
        F.update({item: temp})
        if eval(item) in ItemPart.keys():
            ItemPart[eval(item)].data_i.append(part)
        else:
            node = item_part(item=item)
            node.data_i.append(part)
            ItemPart.update({eval(item): node})


def i_itemset(data, part_minsup, f_1, F):
    for item, t in list(data.items()):
        i_sup = len(t)
        if i_sup < part_minsup:
            del data[item]
        else:
            f_1.append([item, t, i_sup])
            temp = Item()
            temp.item = item
            temp.transet = t
            temp.sup = i_sup
            F.update({item: temp})


def k_Mine(itemlist, flist, F, part_minsup, part):
    st = time.time()
    if '0' in itemlist.keys():  # itemlist={'0':set(flist)}
        item_candidate = {}
        candidate = flist
        for i in range(1, len(candidate)):
            temp_c = set(candidate[:i])
            item_candidate.update({candidate[i]: temp_c})
        k_Mine(item_candidate, flist, F, part_minsup, part)
    else:
        itemset_candidate = {}
        out_candidate = {x: set() for x in flist}
        for i in itemlist.keys():
            item_candidate = {}
            out_candidate.update({i: set()})
            candidate_i = list(itemlist[i])
            candidate_i = sorted(candidate_i, key=flist.index)
            for j in range(len(candidate_i)):
                itemset = i + ',' + candidate_i[j]
                trans = list(set(F[i].transet) & set(F[candidate_i[j]].transet))
                sup = len(trans)
                if sup >= part_minsup:
                    temp = Item()
                    temp.item = itemset
                    temp.transet = trans
                    temp.sup = sup
                    F.update({itemset: temp})
                    temp_c = set(candidate_i[:j]) - out_candidate[i] - out_candidate[candidate_i[j]]
                    if eval(itemset) in ItemPart.keys():
                        ItemPart[eval(itemset)].data_i.append(part)
                    else:
                        node = item_part(item=itemset)
                        node.data_i.append(part)
                        ItemPart.update({eval(itemset): node})
                    if len(temp_c) > 0:
                        item_candidate[itemset] = temp_c
                else:
                    out_candidate[i].add(candidate_i[j])
            itemset_candidate.update(item_candidate)
        if len(itemset_candidate) > 1:
            k_Mine(itemset_candidate, flist, F, part_minsup, part)


def transet(item, part):
    if isinstance(item, int):
        s1 = data[part][str(item)]
    else:
        s1 = data[part][str(item[0])]
        for i in range(1, len(item)):
            s1 = s1 & data[part][str(item[i])]
    return len(s1)


def set_str(s):
    if isinstance(s, int):
        ss = str(s)
    else:
        ss = ''
        for i in s:
            ss = str(i) + ',' + ss
        ss = ss[:-1]
    return ss

if __name__ == '__main__':
    dataset = {
        # 'Chess': [3057, 3017, 2996],
        # 'Connect': [66727, 66423, 66226],
        # 'Kosarak': [32127, 21373, 16389],
        # 'Mushroom': [5040, 4640, 4184],
        # 'Pumsb_star': [32710, 30898, 30204]
        'Chess': [2996]
    }

    # eps_list = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    # topk = [50, 100, 150]
    eps_list = [1.0]
    topk = [150]
    part = {
        'Chess': 2000,
        'Connect': 20000,
        'Kosarak': 100000,
        'Mushroom': 4000,
        'Pumsb_star': 25000
    }
    for database in dataset.keys():
        partlen = part[database]
        for eps in eps_list:
            for sup_i in range(1):
                for j in range(1):
                    infile = 'dataset/' + database + ".txt"
                    t1 = time.time()
                    # Scanning database
                    data, trans, times, last_data_count, MaxTrans = data_read(infile, partlen)
                    t2 = time.time()
                    # print('data reading time use:', t2 - t1)
                    # minsup caculate
                    if dataset[database][sup_i] < 1:
                        minsup = math.floor(dataset[database][sup_i] * trans)
                    else:
                        minsup = dataset[database][sup_i]
                    part_minsup = [math.floor(minsup / trans * partlen)] * times
                    part_minsup.append(math.floor(minsup / trans * last_data_count))
                    supNode = []
                    for i in range(times + 1):
                        supNode.append(data_sup(i))
                    # itemset mining
                    F = {}  # item:sup
                    ItemPart = {}
                    for i in supNode:
                        Fi = partMine(i)
                        F = dict(Counter(F) + Counter(Fi))
                    tm = time.time()
                    # print('mining time use', tm - t2)
                    # itemset sup judge
                    # 隐私预算：privacyBudget=sum(privacyBudget),
                    privacyBudget = eps
                    itemset_count = len(F)
                    lamuda = itemset_count / privacyBudget  # 判断的敏感度

                    fre_itemset = {}
                    for k, v in F.items():
                        n_v = np.random.laplace(0, lamuda)
                        part = set(range(times))
                        part = part - set(ItemPart[k].data_i)
                        sup_k = v
                        if len(part) > 0:
                            if isinstance(k, int):
                                for p in part:
                                    if p in ItemPart[k].data_i:
                                        sup_k = sup_k + transet(k, p)
                            else:
                                item_k = list(k)
                                for p in part:
                                    if all(p in ItemPart[x].data_i for x in item_k):
                                        sup_k = sup_k + transet(item_k, p)
                        n_v = sup_k + n_v
                        # if n_v >= minsup:
                        #     fre_itemset.update({k: n_v})
                        fre_itemset.update({k: n_v})
                    tp = time.time()
                    print('fre_itemset count:', len(fre_itemset))
                    # print('judge time use', tp - tm)
                    # print('total time use:', tp - t1)
                    # 写入outfile
                    fre_item = []
                    for key in fre_itemset.keys():
                        fre_item.append([key, fre_itemset[key]])
                    fre_item_sup = sorted(fre_item, key=lambda x: x[1], reverse=True)
                    s = topk[sup_i]
                    # outfile = '加噪输出\\DP_PartFIM\\' + database + '_eps_' + str(eps) + '_top_' + str(
                    #     s) + '-' + str(j) + '.txt'
                    outfile = database + '_eps_' + str(eps) + '_top_' + str(s) + '-' + str(j) + '.txt'
                    fre_item_k = fre_item_sup[:int(s)]
                    empty_file(outfile)
                    write_into_file('data reading time use:' + str(t2 - t1), outfile)
                    write_into_file('total time use:' + str(tp - t1), outfile)
                    write_into_file('mining file:' + infile, outfile)
                    [write_into_file(set_str(out[0]) + ':' + str(out[1]), outfile) for out in fre_item_k]
                    print("mining time:", tp - t1)
                print(database + " finish mining!")
