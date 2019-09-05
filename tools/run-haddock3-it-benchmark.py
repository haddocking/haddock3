import glob
import os
import time

target_list = glob.glob('*')

n = 2
chunks = [target_list[i:i + n] for i in range(0, len(target_list), n)]

for folder_list in chunks:
	completion_counter = 0
	for f in folder_list:
		if not os.path.isfile(f'{f}/DONE'):
			# submit
			os.system(f'qsub {f}/run.job')
		else:
			completion_counter += 1

	while completion_counter != n:
		# both submitted, wait
		completion_counter = 0
		for f in folder_list:
			if os.path.isfile(f'{f}/DONE'):
				completion_counter += 1
			else:
				time.sleep(60)
