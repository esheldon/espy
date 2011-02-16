#!/usr/bin/env python

from xml.dom import minidom
import sys

def CreateRecipeXML(file):
    """
    Create an example XML recipe doc and see if we can import it
    """

    doc = minidom.Document()

    # Root element
    recipes = doc.createElement("recipes")
    doc.appendChild(recipes)

    # Add a recipe
    recipe = doc.createElement("recipe")

    # Add this node
    recipes.appendChild(recipe)

    # ingredients
    ing = doc.createElement("ingredients")
    recipe.appendChild(ing)

    # Add an item to the recipe
    item = doc.createElement("item")
    ing.appendChild(item)

    # Add the data
    itemdata = doc.createTextNode("Sugar")
    item.appendChild(itemdata)


    # Add an item to the recipe
    item = doc.createElement("item")
    ing.appendChild(item)

    # Add the data
    itemdata = doc.createTextNode("Eggs")
    item.appendChild(itemdata)

    print doc.toprettyxml(indent="    ")

if __name__=="__main__":
    """
    Create an xml document in memory. 
    """

    file = sys.argv[1]
    CreateRecipeXML(file)

