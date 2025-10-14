from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker
from model.base import Write, write_engine, read_engine

if __name__ == "__main__":
    Write.metadata.create_all(write_engine)

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
