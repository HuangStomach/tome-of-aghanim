from typing import List
from typing import Optional

from sqlalchemy import ForeignKey, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Write

class UserCourse(Write):
    __tablename__ = '_r_user_course'
    user_id: Mapped[int] = mapped_column(ForeignKey("user.id"), primary_key=True)
    course_id: Mapped[int] = mapped_column(ForeignKey("course.id"), primary_key=True)
    special_key: Mapped[Optional[str]] = mapped_column(String(50))

    user = relationship("User", back_populates="course_assoc")
    course = relationship("Course", back_populates="user_assoc")

