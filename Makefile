py = python3
clean:
	rm -rf *~ unittest/*~

test:
	$(py) unittest/test_parallel_instability.py