from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import DeclarativeBase

write_meta = MetaData()
read_meta = MetaData()

write_engine = create_engine("mariadb+pymysql://testA:balabala@localhost:3306/testA?charset=utf8mb4", echo=True)
read_engine = create_engine("mariadb+pymysql://testB:balabala@localhost:3307/testB?charset=utf8mb4", echo=True)
# write_engine = create_engine("postgresql+psycopg2://db1", echo=True)
# read_engine = create_engine("postgresql+psycopg2://db2", echo=True)

class Write(DeclarativeBase):
    metadata = write_meta

class Read(DeclarativeBase):
    metadata = read_meta

