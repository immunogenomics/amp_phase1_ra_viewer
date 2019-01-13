$(document).ready(function() {
  document.getElementById("button-data").onclick = function() {
    var url = document.querySelector('a[data-value="data"]').href;
    console.log(url);
    var x = "#" + url.split("#")[1];
    var s = 'a[href="' + x + '"]';
    console.log(s);
    document.querySelector(s).click();
  };
});


