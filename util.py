import os.path
import itertools

# -------------------------------------
#  general util funcs
# -------------------------------------

def close_files(files):
	for file_ in files.itervalues():
		file_.close()

def rewind_files(files):
	for file_ in files.itervalues():
		file_.seek(0)

# return a sorted list of an iterator of elements
def sorted_list(elems):
	l = list(elems)
	l.sort()
	return l

# given a list of unique elements, return a dict of each element and its index in the list
def indices_dict(list_):
	return dict(zip(list_, itertools.count(0)))

# -------------------------------------
#  sample util funcs
# -------------------------------------

def check_files_in_samples(samples):
	all_files_found = True
	for (sample, replicates) in samples.iteritems():
		for replicate in replicates:
			if not os.path.exists(replicate):
				print 'Could not find: %s' % replicate
				all_files_found = False
	return all_files_found

def each_replicate(samples):
	return itertools.chain(*samples.itervalues())
