'''
*----------------------------------*
 Q2DTor and TheRa Programs

 Copyright (c) 2018 Universidade de Santiago de Compostela

 This file is part of both Q2DTor and TheRa softwares.

 Q2DTor and TheRa are free softwares: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Q2DTor and TheRa are distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 inside both Q2DTor and TheRa manuals.  If not, see <http://www.gnu.org/licenses/>.
*----------------------------------*

 This module contains different quotations

*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*
'''

import random, time

def random_quote(program_name):
    quotes = []

    quote  = ""
    quote  = quote + ''' "Your theory is crazy, but it's not crazy enough to be true."'''
    author = "Niels Bohr"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Intelligence is the ability to adapt to change."'''
    author = "Stephen Hawking"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Success is a science; if you have the conditions, you get the result."'''
    author = "Oscar Wilde"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Everything is theoretically impossible, until it is done."'''
    author = "Robert A. Heinlein"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Science is not a body of knowledge, nor a belief system,  \n'''
    quote  = quote + '''  it is just a term which describes human kind's incremental\n'''
    quote  = quote + '''  acquisition of understanding through observation."          '''
    author = "Tim Minchin"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Science is simply the word we use to describe\n'''
    quote  = quote + '''  a method of organizing our curiosity."         '''
    author = "Tim Minchin"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "It took less than an hour to make the atoms,              \n'''
    quote  = quote + '''  a few hundred million years to make the stars and planets,\n'''
    quote  = quote + '''  but five billion years to make man!"                        '''
    author = "George Gamow"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "In some sort of crude sense, which no vulgarity, no humor,           \n'''
    quote  = quote + '''  no overstatement can quite extinguish, the physicists have known sin;\n'''
    quote  = quote + '''  and this is a knowledge which they cannot lose."                       '''
    author = "J. Rober Oppenheimer"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "Why, sir, there is every probability that you will soon be able to tax it."'''
    author = "Michael Faraday"
    quotes.append( (quote,author) )

    quote  = ""
    quote  = quote + ''' "There is nothing that living things do that cannot be  \n'''
    quote  = quote + '''  understood from the point of view that they are made of\n'''
    quote  = quote + '''  atoms acting according to the laws of physics."          '''
    author = "Richard P. Feynman"
    quotes.append( (quote,author) )

    quote, author = random.choice(quotes)
    date = time.strftime("%c")
    return '\n\n%s\n  %s\n\n  End of %s output! Current date & time: %s\n\n'%(quote,author,program_name,date)
