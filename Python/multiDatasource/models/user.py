from typing import List

from sqlalchemy import String, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.associationproxy import AssociationProxy

from models.base import Write

class User(Write):
    __tablename__ = 'user'
    
    id: Mapped[int] = mapped_column('id', primary_key=True)
    name: Mapped[str] = mapped_column('name', String(30))
    email: Mapped[str] = mapped_column('email', String(50))
    age: Mapped[int] = mapped_column('age', nullable=True, index=True)
    gender: Mapped[int] = mapped_column(index=True)

    course_assoc: Mapped[List['UserCourse']] = relationship(back_populates="user")
    courses: AssociationProxy[List['Course']] = association_proxy("course_assoc", 'course')
    
    def __repr__(self):
        return f"<User(name='{self.name}', email='{self.email}')>"

