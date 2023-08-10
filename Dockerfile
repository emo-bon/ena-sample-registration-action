FROM python:3.9
COPY action.py /action.py
COPY entrypoint.sh /entrypoint.sh
COPY requirements.txt /requirements.txt
RUN python -m pip install -r /requirements.txt
ENTRYPOINT ["/entrypoint.sh"]
