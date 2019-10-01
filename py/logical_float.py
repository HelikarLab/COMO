testexp = '(( 2977 ) and ( 2983 )) or 4882 or (( 2982 ) and ( 2974 )) or 2986 or 2984 or 4881 or (( 2977 ) and ( 2974 )) or 3000 or (( 2982 ) and ( 2983 ))'

def float_logical_exp(expressionIn, level=0):
    try:
        loc_r = expressionIn.index(')')
    except:
        if 'and' in expressionIn:
            expressionIn = expressionIn.replace('and', ',')
            expressionIn = 'min{' + expressionIn + '}'
        elif 'or' in expressionIn:
            expressionIn = expressionIn.replace('or', ',')
            expressionIn = 'min{' + expressionIn + '}'
        else:
            expressionIn = expressionIn.replace('[', '')
            expressionIn = expressionIn.replace(']', '')
        return expressionIn
    loc_l = expressionIn[0:loc_r].rindex('(')
    innerstring = expressionIn[loc_l:loc_r+1]
    innerstring = innerstring.replace('(', '[')
    innerstring = innerstring.replace(')', ']')
    if 'and' in innerstring:
        innerstring = innerstring.replace('and', ',')
        innerstring = 'min{'+innerstring+'}'
    elif 'or' in innerstring:
        innerstring = innerstring.replace('or', ',')
        innerstring = 'min{'+innerstring+'}'
    else:
        innerstring = innerstring.replace('[', '')
        innerstring = innerstring.replace(']', '')

    expressionOut = '{}{}{}'.format(expressionIn[0:loc_l], innerstring, expressionIn[loc_r+1:])
    expressionOut = float_logical_exp(expressionOut, level+1)

    return expressionOut

if not (testexp[0] == '(' and testexp[-1] == ')'):
    testexp = '('+testexp+')'
outexp=float_logical_exp(testexp)

outexp = outexp.replace('{','(')
outexp = outexp.replace('}',')')
print(eval(outexp))