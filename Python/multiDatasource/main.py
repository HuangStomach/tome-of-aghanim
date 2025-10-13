from sqlalchemy import create_engine
from sqlalchemy.orm import DeclarativeBase, Session, Mapped, mapped_column

write_meta = MetaData(schema="some_schema")
read_meta = MetaData(schema="some_schema")

class Write(DeclarativeBase):
    metadata = write_meta

class Read(DeclarativeBase):
    metadata = read_meta

class User(Write):
    __tablename__ = 'users'
    
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str]= mapped_column(String(30))
    email = Mapped[str]= mapped_column(String(50))
    
    def __repr__(self):
        return f"<User(name='{self.name}', email='{self.email}')>"

class Book(Read):
    __tablename__ = 'books'
    
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str]= mapped_column(String(30))
    isbn = Mapped[str]= mapped_column(String(50))
    
    def __repr__(self):
        return f"<Book(name='{self.name}', isbn='{self.isbn}')>"

write_engine = create_engine("postgresql+psycopg2://db1", echo=True)
read_engine = create_engine("postgresql+psycopg2://db2", echo=True)
write_meta.create_all(write_engine)

# 创建Session工厂
Session = sessionmaker()
Session.configure(binds={Write: write_engine, Read: read_engine})
session = Session()

# 插入新用户
new_user = User(name='Alice', email='alice@example.com')
session.add(new_user)

# 查询并更新用户
user = session.query(User).filter_by(name='Alice').first()
if user:
    user.email = 'new_alice@example.com'

# 提交事务
session.commit()

# 查询所有用户
users = session.query(User).all()
print(users)

# 关闭会话
session.close()
