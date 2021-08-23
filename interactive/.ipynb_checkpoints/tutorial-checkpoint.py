{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "51275955",
   "metadata": {},
   "source": [
    "# Importing the survey package:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "8700cc87",
   "metadata": {},
   "outputs": [],
   "source": [
    "import sys\n",
    "sys.path.append('..')\n",
    "\n",
    "from survey.protein_builder import get_protein"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "9a72b25a",
   "metadata": {},
   "source": [
    "# Creating a Protein instance for 1Z5F:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "5b6c19b8",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "<survey.proteins.Protein object at 0x7fd8893f01f0>\n"
     ]
    }
   ],
   "source": [
    "pdb_id = '1Z5F'\n",
    "bmrb_id = '5787'\n",
    "\n",
    "protein = get_protein(pdb_id, bmrb_id)\n",
    "print(protein)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "6711bf45",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
