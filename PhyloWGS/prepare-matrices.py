#!/usr/bin/env python2
import argparse
import os
import sys
import cPickle as pickle
import itertools
import numpy as np
import math

def nCk(n, k):
  f = math.factorial
  return f(n) / (f(n - k) * f(k))

def get_trees(tree_dir, num_trees=None):
  fnames = os.listdir(tree_dir)
  fnames = [t for t in fnames if '.DS_Store' not in t]
  fnames = sorted(fnames, key=float)
  fnames = list(reversed(fnames))
  if num_trees is not None:
    fnames = fnames[:num_trees]

  for fname in fnames:
    fpath = os.path.join(tree_dir, fname)
    with open(fpath) as fh:
      tssb = pickle.load(fh)
      score = float(fname)
      yield (score, tssb)

def determine_vertex_relations(tree):
  relations = {}
  wts, all_vertices = tree.get_mixture()

  def _traverse_r(vertex, ancestors):
    for ancestor in ancestors:
      relations[(ancestor, vertex)] = 'ancestor'
      relations[(vertex, ancestor)] = 'descendant'
    for child in vertex.children():
      _traverse_r(child, ancestors + [vertex])

  for vert1, vert2 in itertools.combinations(all_vertices, 2):
    if (vert1, vert2) in relations:
      continue
    relations[(vert1, vert2)] = 'diff_branches'
    relations[(vert2, vert1)] = 'diff_branches'

  _traverse_r(tree.root['node'], [])
  return relations

def generate_ssm_pairs(tree):
  wts, nodes = tree.get_mixture()
  all_ssms = []
  for vertex in nodes:
    ssms = vertex.get_data()
    for ssm in vertex.get_data():
      all_ssms.append((ssm, vertex))
  # Ensure consistent ordering in matrix across different runs.
  all_ssms.sort(key = lambda pair: float(pair[0].name[1:]))
  combos = itertools.combinations(all_ssms, 2)
  num_combos = nCk(len(all_ssms), 2)
  return (num_combos, combos)

def write_ssm_names(tree):
  num_combos, combos = generate_ssm_pairs(tree)
  with open('ssm_pair_names', 'w') as ssm_names:
    for (ssm1, vertex1), (ssm2, vertex2) in combos:
      ssm_names.write('%s,%s\n' % (ssm1.name, ssm2.name))

def create_matrix(tree):
  num_combos, combos = generate_ssm_pairs(tree)
  matrix = np.zeros((4, num_combos), dtype=np.int8)
  vertex_relations = determine_vertex_relations(tree)

  i = 0
  for (ssm1, vertex1), (ssm2, vertex2) in combos:
    if vertex1 == vertex2:
      # SSMs occur in same population
      ssm_relation = 0
    else:
      vert_relation = vertex_relations[vertex1, vertex2]
      if vert_relation == 'ancestor':
        ssm_relation = 1
      elif vert_relation == 'descendant':
        ssm_relation = 2
      elif vert_relation == 'diff_branches':
        ssm_relation = 3
      else:
        raise Exception('Unknown SSM relation between %s and %s: %s' % (vert1, vert2, vert_relation))
    matrix[ssm_relation, i] = 1
    i += 1

  return matrix

def calculate_matrices(tree_dir):
  idx = 0
  first_iteration = True

  for score, tree in get_trees(tree_dir):
    matrix = create_matrix(tree)
    if first_iteration:
      mean_matrix = np.zeros(matrix.shape, dtype=np.float64)
      # Only need to write these on the first iteration, as they will be
      # identical for all trees.
      write_ssm_names(tree)
      first_iteration = False

    mean_matrix += matrix
    outfn = 'mutpairs_%s.txt.gz' % score
    #np.savetxt(outfn, matrix)

    idx += 1

  mean_matrix /= idx
  np.savetxt('mutpairs_mean.txt.gz', mean_matrix)

def main():
  parser = argparse.ArgumentParser(description='Plot stats about each tree')
  parser.add_argument('tree_dir')
  args = parser.parse_args()
  calculate_matrices(args.tree_dir)

main()
