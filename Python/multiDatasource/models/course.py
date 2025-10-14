from typing import List

from sqlalchemy import String, Integer
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.associationproxy import AssociationProxy

from models.base import Write

class Course(Write):
    __tablename__ = 'course'
    
    id: Mapped[int] = mapped_column('id', primary_key=True)
    name: Mapped[str] = mapped_column('name', String(30))

    user_assoc: Mapped[List['UserCourse']] = relationship(back_populates="course")
    users: AssociationProxy[List['User']] = association_proxy("user_assoc", "user")
    
    def __repr__(self):
        return f"<Course(name='{self.name}')>"

