.PHONY: test lint install clean all run demo full_analysis

full_analysis:
	python run_full_analysis.py

demo:
	python demo.py

run:
	python app.py --port 5001

all: install lint test

clean:
	find . -type f -name '*.pyc' -delete
	find . -type d -name '__pycache__' -delete

install:
	pip install -r requirements.txt

test:
	python -m unittest discover tests

lint:
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics 