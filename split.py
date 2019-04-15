# -*- coding:utf-8 -*-
from datetime import datetime


def Main():
    # source_dir = '/mnt/hgfs/transopon_data/data/output/bowtie_l1hs/A152-Normal.s5_9_Sl31275..fastq.cleaned.fastq.sam'
    # target_dir = '/mnt/hgfs/transopon_data/data/output/bowtie_l1hs/split/'
    source_dir ='/mnt/hgfs/transopon_data/data/output/bowtie_human/A152-Normal.s5_9_Sl31275..fastq.cleaned.fastq.sam'
    target_dir = '/mnt/hgfs/transopon_data/data/output/bowtie_human/'
    # 计数器
    flag = 0

    # 文件名
    name = 1

    # 存放数据
    datalist = []

    print("开始。。。。。")
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    with open(source_dir, 'r') as f_source:
        for line in f_source:
            flag += 1
            datalist.append(line)
            if flag == 50000:
                with open(target_dir + "sam" + str(name) + ".txt", 'w+') as f_target:
                    for data in datalist:
                        f_target.write(data)
                name += 1
                flag = 0
                datalist = []
                break
    # 处理最后一批行数少于200万行的
    # with open(target_dir + "sam" + str(name) + ".txt", 'w+') as f_target:
    #     for data in datalist:
    #         f_target.write(data)

    print("完成。。。。。")
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == "__main__":
    Main()
