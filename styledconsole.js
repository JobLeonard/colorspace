// Based on http://jsfiddle.net/yg6hk/5/

const styledConsoleLog = (() => {
  // Initiate regex once and re-use
  var startTagRe = /\{\{(['"])([^'"]*)\1\s*/gi;
  var endTagRe = /\}\}/gi;
  return function () {
    var argArray = [];

    if (arguments.length) {

      var reResultArray;
      argArray.push(arguments[0].replace(startTagRe, '%c').replace(endTagRe, '%c'));
      while (reResultArray = startTagRe.exec(arguments[0])) {
        argArray.push(reResultArray[2]);
        argArray.push('');
      }

      // pass through subsequent args since chrome dev tools does not (yet)
      // support console.log styling of the following form:
      // console.log('%cBlue!', 'color: blue;', '%cRed!', 'color: red;');
      for (var j = 1; j < arguments.length; j++) {
        argArray.push(arguments[j]);
      }
    }

    console.log.apply(console, argArray);
  };
})();

styledConsoleLog(`
{{"color:hsl(0, 100%, 90%);background-color:hsl(0, 100%, 50%);" Red }}
{{"color:hsl(39, 100%, 85%);background-color:hsl(39, 100%, 50%);" Orange }}
{{"color:hsl(60, 100%, 35%);background-color:hsl(60, 100%, 50%);" Yellow }}
{{"color:hsl(120, 100%, 60%);background-color:hsl(120, 100%, 25%);" Green }}
{{"color:hsl(240, 100%, 90%);background-color:hsl(240, 100%, 50%);" Blue }}
{{"color:hsl(300, 100%, 85%);background-color:hsl(300, 100%, 25%);" Purple }}
{{"color:hsl(0, 0%, 80%);background-color:hsl(0, 0%, 0%);" Black }}
`);


var printColors = (() => {
  let blockString = '█████████████████████████████████';
  let pc =  (colorArray) => {
    let colorCSS = [],
      colorString = [],
      idx1 = colorArray.length;
    while (idx1--) {
      let idx2 = colorArray.length - 1 - idx1,
        idx3 = (idx1 + (colorArray.length >> 1)) % colorArray.length,
        idx4 = (idx2 + (colorArray.length >> 1)) % colorArray.length,
        hex1 = colorArray[idx1],
        hex2 = colorArray[idx2],
        hex3 = colorArray[idx3],
        hex4 = colorArray[idx4];

      let idxString = idx1.toString(10),
        idxString2 = idx2.toString(10),
        idxString3 = idx3.toString(10),
        idxString4 = idx4.toString(10);
      while (idxString.length < 3) {
        idxString = '0' + idxString;
      }
      while (idxString2.length < 3) {
        idxString2 = '0' + idxString2;
      }
      while (idxString3.length < 3) {
        idxString3 = '0' + idxString3;
      }
      while (idxString4.length < 3) {
        idxString4 = '0' + idxString4;
      }
      colorCSS.push(
        'font-weight: bold;',
        `color:${hex1}; font-weight: bold;`,
        'font-weight: bold;',
        `color:${hex2}; font-weight: bold;`,
        'font-weight: bold;',
        `color:${hex3}; font-weight: bold;`,
        'font-weight: bold;',
        `color:${hex4}; font-weight: bold;`,

        `color:${hex1}; font-weight: bold;`,
        `color:${hex2}; font-weight: bold;`,
        `color:${hex3}; font-weight: bold;`,
        `color:${hex4}; font-weight: bold;`,
      );
      colorString.push(
        `
%c${idxString}:%c${blockString}%c ${idxString2}:%c${blockString}%c ${idxString3}:%c${blockString}%c ${idxString4}:%c${blockString}`,`
%c    ${blockString}%c     ${blockString}%c     ${blockString}%c     ${blockString}`
      );
    }
    console.log(colorString.join(''), ...colorCSS);
  };

  // test
  pc(['#ffffff', '#ff0000', '#00ff00', '#0000ff']);
  return pc;
})();



function swap(arr, i, j) {
  arr = arr.slice(0);
  let t = arr[i];
  arr[i] = arr[j];
  arr[j] = t;
  printColors(arr);
  return arr;
}