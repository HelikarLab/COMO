from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
#Let’s import our Book and Base classes from our database_setup.py file
from database_setup import Book, Base

engine = create_engine('sqlite:///books-collection.db')
# Bind the engine to the metadata of the Base class so that the
# declaratives can be accessed through a DBSession instance
Base.metadata.bind = engine

DBSession = sessionmaker(bind=engine)
# A DBSession() instance establishes all conversations with the database
# and represents a "staging zone" for all the objects loaded into the
# database session object.
session = DBSession()

