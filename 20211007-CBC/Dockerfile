FROM python:3.8.2

ADD ./2_submit2_self_addid.py /code/2_submit2_self_addid.py
ADD ./requirements.txt /code/requirements.txt
ADD ./trained_model.pth /code/trained_model.pth

WORKDIR /code

RUN pip install -r requirements.txt -i https://mirrors.aliyun.com/pypi/simple/
