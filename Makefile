install:
	pip install -e .[dev]

lint:
	flake8 just_psf tests --max-line-length=120 --ignore=N802

test:
	pytest tests