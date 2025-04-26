test:
	pytest -r a -v tests --cov mzspeclib --cov-report=html --cov-report term

retest:
	py.test -v tests --lf --pdb


rebuild_test_bundle:
	mzspeclib convert -f text tests/test_data/chinese_hamster_hcd_selected_head.msp tests/test_data/chinese_hamster_hcd_selected_head.mzspeclib.txt
	mzspeclib convert -f json tests/test_data/chinese_hamster_hcd_selected_head.msp tests/test_data/chinese_hamster_hcd_selected_head.mzspeclib.json
	mzspeclib convert -f json "tests/test_data/complex_interpretations_with_members.mzspeclib.txt" tests/test_data/complex_interpretations_with_members.mzspeclib.json
	python tests/test_data/generate_annotations.py
	mzspeclib convert -f text examples/chinese_hamster_hcd_selected_head.msp examples/chinese_hamster_hcd_selected_head.mzspeclib.txt

update-readme:
	cog -r -U ./README.md