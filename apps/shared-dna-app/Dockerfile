FROM python:3.10-slim

USER root

WORKDIR /app

ADD . /app

RUN pip install --trusted-host pypi.python.org -r requirements.txt

EXPOSE 8001

ENV NAME shared-dna

CMD ["gunicorn", "-b", "0.0.0.0:8001", "-t", "600", "app:app"]
