
QUnit.module('RecipesTest', {
    beforeEach: function() {
        //var $ = django.jQuery;
    }
});

QUnit.test( "Test basic AJAX functionality in recipes", function( assert ) {
  //assert.ok( true, "Passed!" );
  //var $ = jQuery;
  //let add_snippet = $('.cmd-value').click();
  assert.dom('.cmd-value').exists();
  //assert.dom(add_snippet.length, 1, 'Error adding snippets to template ')

});