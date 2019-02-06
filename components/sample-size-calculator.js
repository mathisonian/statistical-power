const React = require('react');
const { doSizeAssoc } = require('./sizelib');

class CustomComponent extends React.Component {

  constructor(props) {
    super(props);
  }


  render() {
    const { hasError, idyll, updateProps, alpha, power, effectSize, ...props } = this.props;
    const calc = doSizeAssoc(alpha, undefined, undefined, undefined, power, effectSize);

    return (
      <span {...props}>{calc.n}</span>
    );
  }
}

module.exports = CustomComponent;
