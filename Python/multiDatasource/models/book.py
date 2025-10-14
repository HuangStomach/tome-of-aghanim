from sqlalchemy import String
from sqlalchemy.orm import Mapped, mapped_column

from models.base import Read

class Book(Read):
    __tablename__ = 'book'
    
    id: Mapped[int] = mapped_column('id', primary_key=True)
    name: Mapped[str] = mapped_column('name', String(30))
    isbn: Mapped[str] = mapped_column('isbn', String(50))
    
    def __repr__(self):
        return f"<Book(name='{self.name}', isbn='{self.isbn}')>"

