#!/usr/bin/env python3
#-*-encoding: utf-8 -*-


def get_w(arq,n):
	W = []
	for i in range(n):
		line = arq.readline()
		W.append([float(x) for x in line.split()[1:]])

	return W


def average(w):
	sum_w = 0
	total_w = 0

	for scene in w:
		sum_w += sum(scene)
		total_w += len(scene)

	return sum_w/total_w


n = 4

file_cg = "calibration/Rect10x10_km/results_cg_%d.txt" % n
file_queue = "calibration/Rect10x10_km/results_queue_%d.txt" % n

arq_cg = open(file_cg, "r")
arq_queue = open(file_queue, "r")

cg = get_w(arq_cg, n)
queue = get_w(arq_queue,n)

avg_cg = average(cg)
avg_queue = average(queue)



print(avg_cg,avg_queue)


arq_cg.close()
arq_queue.close()