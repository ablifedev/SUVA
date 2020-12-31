#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

####################################################################################
# Copyright (C) 2015-2019 by ABLIFE
####################################################################################
# 名称：expression_quantity_calculation.py
# 描述：计算表达量
# 作者：程超
# 创建时间：2015-4-10
# 联系方式：chaocheng@ablife.cc
####################################################################################
# 修改记录
####################################################################################
# Date           Version       Author            ChangeLog
# 2015-4-10      v0.1          ChengChao         创建测试版本
# 2015-7-24      v1.0          ChengChao         加入对bed文件支持；对于bed格式先扫描所有
#                                                所有基因，记录区域信息，再扫描bed iv
#
#
#####################################################################################

"""
程序功能说明：
1.计算gene表达量
2.randCheck_gene
3.randCheck_mRNA
程序设计思路：
利用gffutils和HTSeq包进行统计
"""

# 导入必要的包
import re
import os
import sys
import logging
import time
import datetime

from optparse import OptionParser, OptionGroup

reload(sys)
sys.setdefaultencoding('utf-8')

import subprocess
import threading
import gffutils
import HTSeq
import numpy
import multiprocessing
import signal
from matplotlib import pyplot

sys.path.insert(1, os.path.split(os.path.realpath(__file__))[0] + "/../")
# print(sys.path)
from ablib.utils.tools import *

# 检查python的版本，我们需要使用python2.7
# TODO: HTSeq升级到python3版本后升级程序到python3
if sys.version_info < (2, 7):
    print("Python Version error: please use phthon2.7")
    sys.exit(-1)

# 程序版本号
_version = 'v0.1'


# -----------------------------------------------------------------------------------
# --- S 参数设置模块
# -----------------------------------------------------------------------------------
def configOpt():
    """Init for option
    """
    usage = 'Usage: %prog [-f] [other option] [-h]'
    p = OptionParser(usage)
    # basic options
    p.add_option('-t', '--totalsj', dest='totalsj',
                 action='store', type='string', help='totalsj file')
    p.add_option('-b', '--brfile', dest='brfile',
                 action='store', type='string', help='br file')
    p.add_option('-o', '--outfile', dest='outfile', default='Mapping_distribution.txt',
                 action='store', type='string', help='gene expression file')

    group = OptionGroup(p, "Preset options")
    # preset options
    group.add_option('-O', '--outDir', dest='outDir', default='./',
                     action='store', type='string', help='output directory', metavar="DIR")
    group.add_option('-L', '--logDir', dest='logDir', default='', action='store',
                     type='string', help='log dir ,default is same as outDir')
    group.add_option('-P', '--logPrefix', dest='logPrefix', default='',
                     action='store', type='string', help='log file prefix')
    group.add_option('-E', '--email', dest='email', default='none', action='store', type='string',
                     help='email address, if you want get a email when this job is finished,default is no email', metavar="EMAIL")
    group.add_option('-Q', '--quiet', dest='quiet', default=False,
                     action='store_true', help='do not print messages to stdout')
    group.add_option('-K', '--keepTemp', dest='keepTemp',
                     default=False, action='store_true', help='keep temp dir')
    group.add_option('-T', '--test', dest='isTest', default=False,
                     action='store_true', help='run this program for test')
    p.add_option_group(group)

    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)

    opt, args = p.parse_args()
    return (p, opt, args)


def listToString(x):
    """获得完整的命令
    """
    rVal = ''
    for a in x:
        rVal += a + ' '
    return rVal


# 解析参数
opt_parser, opt, args = configOpt()

if opt.logDir == "":
    opt.logDir = opt.outDir + '/log/'


# 对参数进行有效性验证和初步处理


# -----------------------------------------------------------------------------------
# --- E 参数设置模块
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# S (全局)变量定义及初始化设置模块
# -----------------------------------------------------------------------------------
# 获取路径信息
scriptPath = os.path.abspath(os.path.dirname(__file__))  # absolute script path
binPath = scriptPath + '/bin'  # absolute bin path
outPath = os.path.abspath(opt.outDir)  # absolute output path
os.mkdir(outPath) if not os.path.isdir(outPath) else None
logPath = os.path.abspath(opt.logDir)
os.mkdir(logPath) if not os.path.isdir(logPath) else None
tempPath = outPath + '/temp/'  # absolute bin path
# os.mkdir(tempPath) if not os.path.isdir(tempPath) else None
resultPath = outPath + '/result/'


# os.mkdir(resultPath) if not os.path.isdir(resultPath) else None

# -----------------------------------------------------------------------------------
# E (全局)变量定义及初始化设置模块
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 日志模块logging初始化
# -----------------------------------------------------------------------------------
def initLogging(logFilename):
    """Init for logging
    """
    logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s : %(levelname)s] %(message)s',
                        datefmt='%y-%m-%d %H:%M', filename=logFilename, filemode='w')
    if not opt.quiet:  # 非quiet模式在屏幕打印出程序执行INFO
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter(
            '[%(asctime)s : %(levelname)s] %(message)s', datefmt='%y-%m-%d %H:%M')
        # tell the handler to use this format
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)


dt = datetime.datetime.now()
logFile = logPath + '/' + opt.logPrefix + 'log.' + \
    str(dt.strftime('%Y%m%d.%H%M%S.%f')) + '.txt'
initLogging(logFile)
logging.debug(sys.modules[__name__].__doc__)  # 打印出程序的说明文档
# -----------------------------------------------------------------------------------
# E 日志模块logging初始化
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 计时器模块->Getting Start Time
# -----------------------------------------------------------------------------------
logging.debug('Program version: %s' % _version)
logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = datetime.datetime.now()
logging.debug("计时器：Program start at %s" % startTime)


# -----------------------------------------------------------------------------------
# E 计时器模块->Getting Start Time
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 类定义（若有）
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# E 类定义（若有）
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 功能函数定义
# -----------------------------------------------------------------------------------
# def invert_strand(iv):
#     """
#     :param iv: HTSeq.GenomicInterval object
#     :return: HTSeq.GenomicInterval - strand is reversed
#     """
#     iv2 = iv.copy()
#     if iv2.strand == "+":
#         iv2.strand = "-"
#     elif iv2.strand == "-":
#         iv2.strand = "+"
#     else:
#         raise ValueError, "Illegal strand"
#     return iv2

def getTotalBase(iv, coverage):
    totalbases = 0
    for iv2, value2 in coverage[iv].steps():
        if value2 > 0:
            totalbases += value2 * iv2.length
    return totalbases


# @profile
def readChrwithBam():
    # print(chr)
    reads_dict = {}

    ss_complex = {}
    ss_complex_total = {}
    ss_use = {}
    

    totalsjfile = opt.totalsj

    # gas = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    ga = HTSeq.GenomicArray("auto", stranded=True, typecode='i')

    f = open(totalsjfile)

    for eachLine in f:
        line = eachLine.strip().split("\t")

        # if int(line[4])<opt.sjreads:
        #     continue
        ch = line[0]
        s = int(line[1])
        e = int(line[2])

        strand = "+"
        if line[3] == "2":
            strand = "-"

        r = int(line[6])

        iv1 = HTSeq.GenomicInterval(ch, s-1, e, strand)

        ga[iv1] += r

        key = ch + "\t" + str(s) + "\t" + strand
        key2 = ch + "\t" + str(e) + "\t" + strand

        if key not in ss_complex:
            ss_complex[key] = list()
        ss_complex[key].append(str(e)+":"+line[6])
        if key2 not in ss_complex:
            ss_complex[key2] = list()
        ss_complex[key2].append(str(s)+":"+line[6])

        if key not in ss_complex_total:
            ss_complex_total[key] = 0
        ss_complex_total[key]+=int(line[6])
        if key2 not in ss_complex_total:
            ss_complex_total[key2] = 0
        ss_complex_total[key2]+=int(line[6])

    f.close()

    f2 = open(totalsjfile)

    for eachLine in f2:
        # print(eachLine)
        line = eachLine.strip().split("\t")

        # if int(line[4])<opt.sjreads:
        #     continue

        ch = line[0]
        s = int(line[1])
        e = int(line[2])

        strand = "+"
        if line[3] == "2":
            strand = "-"

        ip = HTSeq.GenomicPosition(ch, s-1, strand)
        ip2 = HTSeq.GenomicPosition(ch, e-1, strand)

        key = ch + "\t" + str(s) + "\t" + strand
        key2 = ch + "\t" + str(e) + "\t" + strand

        reads_dict[key]=ga[ip]
        reads_dict[key2]=ga[ip2]

        

        if key not in ss_use:
            ss_use[key] = list()
        ratio = round(float(line[6])/ss_complex_total[key2],2)
        ss_use[key].append(str(e)+":"+str(ratio))
        if key2 not in ss_use:
            ss_use[key2] = list()
        ratio = round(float(line[6])/ss_complex_total[key],2)
        ss_use[key2].append(str(s)+":"+str(ratio))

    f2.close()
    return reads_dict, ss_complex, ss_use
    

    # del reads_dict

    # logging.info("done %s" % chr)


# -----------------------------------------------------------------------------------
# E 功能函数定义
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 主函数
# -----------------------------------------------------------------------------------
def main():
    print("Main procedure start...")

    # 1.读入gff文件/gtf文件/annotationDB
    # 读取gff，建立database，以info中的ID信息作为gene的标识，如果db文件已经存在则不需要再提供gff选项，否则db文件会被覆盖
    # chrs = {}
    # for chr in os.popen("cut -f 1 " + bedfile + " | sort |uniq").readlines():
    #     chr = chr.strip()
    #     # print(chr)
    #     chrs[chr] = 1

    # 2.对每个染色体多线程处理，遍历每个gene，读取gene内的reads，进行计算

    # for chr in chrs:
    #     if chr == "chrM":
    #         continue
    #     if not chr.startswith("chr"):
    #         continue
    #     reads[chr] = {}
    #     readChrwithBam(chr, reads)
    
    reads_dict,ss_complex,ss_use = readChrwithBam()

    br_dict = {}
    f = open(opt.brfile)
    for eachLine in f:
        # chr10   100042573   100048757   -   172 1   0
        line = eachLine.strip().split("\t")

        ch = line[0]
        s = int(line[1])+1
        e = int(line[2])

        strand = line[3]
        key = ch + "\t" + str(s) + "\t" + strand
        key2 = ch + "\t" + str(e) + "\t" + strand
        br_dict[key] = int(line[5])
        br_dict[key2] = int(line[6])
    f.close()

    w = open(opt.outfile, 'w')
    w.writelines("chr\tpos\tstrand\tsjcoverage\tfrom:sj\tfrom_num\tto:ratio\n")

    for s in sorted(reads_dict):
        if s not in br_dict:
            br_dict[s]=0
        ssc = ",".join(ss_complex[s])
        sslen = len(ss_complex[s])
        suc = ",".join(ss_use[s])
        # sulen = len(ss_use[s])
        w.writelines(s + "\t"+str(reads_dict[s]+br_dict[s])+"\t"+ssc+"\t"+str(sslen)+"\t"+suc+"\n")

    w.close()


if __name__ == '__main__':
    # 执行主函数
    main()
# -----------------------------------------------------------------------------------
# E 主函数
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 清理临时文件夹，放在最后
# -----------------------------------------------------------------------------------
if not opt.keepTemp:
    os.system('rm -rf ' + tempPath)
    logging.debug("Temp folder is deleted..")
# -----------------------------------------------------------------------------------
# E 清理临时文件夹，放在最后
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 计时器模块->Getting Total Run Time
# -----------------------------------------------------------------------------------
logging.debug("Program ended")
currentTime = datetime.datetime.now()
runningTime = (currentTime - startTime).seconds  # in seconds
logging.debug("计时器：Program start at %s" % startTime)
logging.debug("计时器：Program end at %s" % currentTime)
logging.debug("计时器：Program ran %.2d:%.2d:%.2d" %
              (runningTime / 3600, (runningTime % 3600) / 60, runningTime % 60))
# -----------------------------------------------------------------------------------
# S 计时器模块->Getting Total Run Time
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# S 发送邮件模块
# -----------------------------------------------------------------------------------
if opt.email != "none":
    run_cmd = listToString(sys.argv)
    sendEmail(opt.email, str(startTime), str(currentTime), run_cmd, outPath)
    logging.info("发送邮件通知到 %s" % opt.email)


# -----------------------------------------------------------------------------------
# S 发送邮件模块
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# S 程序运行计数器
# -----------------------------------------------------------------------------------
def countProgram(programName, startT, runT, isTest):
    countProgramFile = open('/users/ablife/ablifepy/countProgram.txt', 'a')
    countProgramFile.write(programName + '\t' + str(os.getlogin()) +
                           '\t' + str(startT) + '\t' + str(runT) + 's\t' + isTest + '\n')
    countProgramFile.close()


testStr = 'P'
if opt.isTest:
    testStr = 'T'
countProgram(sys.argv[0], startTime, runningTime, testStr)
# -----------------------------------------------------------------------------------
# E 程序运行计数器
# -----------------------------------------------------------------------------------
