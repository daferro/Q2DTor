'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Q2DTor
Version     : 1.1
License     : MIT/x11

Copyright (c) 2019, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------


*----------------------------------*
| Main Author:  David Ferro-Costas |
| Last Update:  Mar-02st-2018      |
*----------------------------------*

 This module contains different quotations

'''

import random, time

def random_quote(program_name):
    quotes = []

    quote  = ''' "Your theory is crazy, but it's not crazy enough to be true."'''
    author = "Niels Bohr"
    quotes.append( (quote,author) )

    quote  = ''' "Intelligence is the ability to adapt to change."'''
    author = "Stephen Hawking"
    quotes.append( (quote,author) )

    quote  = ''' "Success is a science; if you have the conditions, you get the result."'''
    author = "Oscar Wilde"
    quotes.append( (quote,author) )

    quote  = ''' "Everything is theoretically impossible, until it is done."'''
    author = "Robert A. Heinlein"
    quotes.append( (quote,author) )

    quote  = ''' "Science is not a body of knowledge, nor a belief system,  \n'''
    quote += '''  it is just a term which describes human kind's incremental\n'''
    quote += '''  acquisition of understanding through observation."          '''
    author = "Tim Minchin"
    quotes.append( (quote,author) )

    quote  = ''' "Science is simply the word we use to describe\n'''
    quote += '''  a method of organizing our curiosity."         '''
    author = "Tim Minchin"
    quotes.append( (quote,author) )

    quote  = ''' "It took less than an hour to make the atoms,              \n'''
    quote += '''  a few hundred million years to make the stars and planets,\n'''
    quote += '''  but five billion years to make man!"                        '''
    author = "George Gamow"
    quotes.append( (quote,author) )

    quote  = ''' "In some sort of crude sense, which no vulgarity, no humor,           \n'''
    quote += '''  no overstatement can quite extinguish, the physicists have known sin;\n'''
    quote += '''  and this is a knowledge which they cannot lose."                       '''
    author = "J. Rober Oppenheimer"
    quotes.append( (quote,author) )

    quote  = ''' "Why, sir, there is every probability that you will soon be able to tax it."'''
    author = "Michael Faraday"
    quotes.append( (quote,author) )

    quote  = ''' "There is nothing that living things do that cannot be  \n'''
    quote += '''  understood from the point of view that they are made of\n'''
    quote += '''  atoms acting according to the laws of physics."          '''
    author = "Richard P. Feynman"
    quotes.append( (quote,author) )

    quote, author = random.choice(quotes)
    date = time.strftime("%c")
    return '\n\n%s\n  %s\n\n  End of %s output! Current date & time: %s\n\n'%(quote,author,program_name,date)
