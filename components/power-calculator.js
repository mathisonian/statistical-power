const React = require('react');
const { doSizeAssoc } = require('./sizelib');

class CustomComponent extends React.Component {

  constructor(props) {
    super(props);
  }


  render() {
    const { hasError, idyll, updateProps, alpha, sampleSize, effectSize, ...props } = this.props;
    const calc = doSizeAssoc(alpha, undefined, undefined, sampleSize, undefined, effectSize);
    return (
      <span {...props} style={{fontFamily:'monospace'}}>{calc.prob}%</span>
    );
  }
}

module.exports = CustomComponent;
