from biostar import VERSION
from django.conf import settings
from biostar.recipes.const import *


def engine(request):
    '''
    Additional context applied to each request.
    Note: This function is critically important!
    The site will not load up without it.
    '''

    params = dict(user=request.user, VERSION=VERSION, request=request,
                  data_board=DATA_BOARD, data_specific=DATA_SPECIFIC_BOARD,
                  recipe_board=RECIPE_BOARD, file_board=FILES_BOARD,
                  recipe_link=LINKED_RECIPES, results_board=RESULTS_BOARD,
                  recipe_specific=RECIPE_SPECIFIC_BOARD,
                  )

    return params
